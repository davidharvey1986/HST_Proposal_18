
'''
New idea to reshape the proposal.

Get error in beta from estaimtes that include
-N halos
-Contamination

This beta is then converted to sigma via the fitted relation in HST constraints

'''
import pickle as pkl
from GetAndPlotTrend import *
from matplotlib import pyplot as plt
from scipy.stats import norm
from HSTconstraints import * #where beta2cross is
def main():
    #contamination fraction sest from Harvey13
    #frac=[0.01,0.1,0.2,0.3]
    #contaminationImprovement = [17./17.,24/17.,30/17.,36/17.]
    #frac=np.linspace(0.01,0.3,10)
    #assume that H15 had a contamunation of 30%

    contaminationImprovement= np.array([0.01,0.25,0.3])
    contaminationLabel= \
      np.array(['Single Band (30% Contamination)',\
                'Dual Band  (10% Contamination)', \
                'Multi Band (1% Contamination)'])
        
    #error change due to the contamintaion improvement measured
    #from Harvey 13
    errorChangeDueToContamination = \
      (contaminationImprovement*60.+17)/17.
    nClusters = np.linspace(20.,300,29)
    ax = plt.gca()
    colors = ['red','orange','green']

    nHalosPerCluster = 2
    for jFrac, jCont in enumerate(errorChangeDueToContamination):
        limits=[]
        print jCont
        for i in nClusters:
            
            iLimit = \
              getCrossSectionLimit(nClusters = i, \
                        contamination=jCont)
            limits.append(iLimit)

        
        ax.plot(nClusters/nHalosPerCluster, \
                    limits, color=colors[jFrac],\
                    label=contaminationLabel[jFrac])

    ax.text(100.,0.2,'Required limits',fontsize=15)
    ax.plot([0.,150],[0.2,0.2],'k--')
    ax.set_ylabel(r'Predicted limit on cross-section cm$^2$/g')
    ax.set_xlabel(r'Required number of merging clusters')
    ax.legend(title='Cluster Galaxy Member \n       Contamination')
    ax.set_xlim(20,140)

    plt.savefig('../plots/requiredNumberClusters.pdf')
    plt.show()

    
def getCrossSectionLimit(   nClusters = 1, contamination=1.):
    #First get the data from the simulations
    crossSection, totBeta, errBeta, totalSI, totalSG = \
      pkl.load(open('Pickles/SimulationBetaToCross.pkl','rb'))
      
    BetaCrossSectionTrend,BetaCrossSectionTrendError = \
      getTrend( crossSection[1:], totBeta[1:], None )

    #I need to assume a beta that is non zero so make it small
    assumedBeta = 0.
    #Have some error that goes as the sqrt of the number of halos
    errorSI = 23./np.sqrt(nClusters)/contamination
    #assume that SG is that from the simulations
    assumedDistSG = np.mean(totalSG[1:])
    #so corresponding beta dist SI is
    assumedDistSI = assumedBeta*assumedDistSG
    #And the error on that beta would be
    errorOnBeta = getErrorOnBeta( errorSI, assumedDistSI, assumedDistSG)
    #Now create a PDF of that beta by first taking an x axis linear in
    #cross section
    nCross=10000
    plotCrossSections = np.linspace(0.,10.,nCross)
    #Beta cannot be greater than Aparam since the log must be finite
    #i.e. beta > BetaCrossSectionTrend[0]
    validBeta = BetaCrossSectionTrend[0]
    plotBetas = np.linspace(0.,validBeta,nCross)[:-1] 
    
    #And create the Gaussian PDF
    pdfBeta = norm.pdf( plotBetas, assumedBeta, errorOnBeta)
    pdfBeta /= np.sum(pdfBeta*plotBetas[1])
    
    #Now convert the pdf into pdf on cross-section
    crossSection = betaToCross( plotBetas, *BetaCrossSectionTrend)
    
    #And times by the jacobian
    derivBetaDerivCross = BetaCrossSectionTrend[1]
    derivBetaDerivCross = \
      getDerivBetaDerivCrossBeta( crossSection, *BetaCrossSectionTrend)

    pdfCross = pdfBeta*derivBetaDerivCross

    #create a cumulative probability and normalise
    dX = crossSection - np.roll(crossSection,1) 
    
    cdfCross = np.cumsum(pdfCross[:-1]*dX[1:])
    
    
    return crossSection[1:][np.argmin(np.abs(cdfCross-0.68))]

def getErrorOnBeta( error_si, delta_si, delta_sg ):
    x = delta_si / delta_sg
    return error_si/delta_sg*np.sqrt(1.+x**2)

def getDerivBetaDerivCrossBeta( x, p1, p2):
    return p1/p2*np.exp(-x/p2)

def testPlotTrend():
    #First get the data from the simulationsa
    crossSection, totBeta, errBeta, totalSI, totalSG = \
      pkl.load(open('Pickles/SimulationBetaToCross.pkl','rb'))
      
    BetaCrossSectionTrend,BetaCrossSectionTrendError = \
      getTrend( crossSection, totBeta, None )


                  
    ax = plt.gca()
    color='blue'
    xPlot = np.linspace(0.,1., 100)
    yPlot = fitfunc( xPlot, *BetaCrossSectionTrend )
    ax.plot( xPlot, yPlot, '--', color=color)

    ax.errorbar( crossSection, totBeta, errBeta)
    trendParamsWithError = np.append( BetaCrossSectionTrend, \
                            BetaCrossSectionTrendError )
                            
    estTrendError = trendError( xPlot, *trendParamsWithError)
    
    ax.fill_between( xPlot, yPlot+estTrendError/2., \
                            yPlot-estTrendError/2., \
                         alpha=0.3, color=color )

    ax.set_xscale('log')
    plt.show()

def trendError( x, p1, p2, ep1, ep2):
    dP1dx = (1.-np.exp(-x/p2))
    dP2dx = p1*x/p2**2*np.exp(-x/p2)
    return np.sqrt( ep1**2*dP1dx**2 + ep2**2*dP2dx**2)
