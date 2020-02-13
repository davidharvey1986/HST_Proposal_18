from getMergerSample import *
import getSimData as GSD
from matplotlib import pyplot as plt

from scipy.stats import norm
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import cauchy
from scipy.stats import rayleigh
from GetAndPlotTrend import *
import matplotlib as mpl

def HSTconstraints():

    '''
    plot and devise the constraints for the HST proposal
    '''
    SimNames = ['CDM','SIDM0.1','SIDM0.3','SIDM1'] #GSD.getSimNameList()
    fig, axarr = plt.subplots( 1 )

    bins=40
    colors = ['k','r','b']
    mpl.rcParams['lines.linewidth'] = 2

    totBeta = np.zeros( len(SimNames))
    stdTotBeta = np.zeros( len(SimNames))
    totDistSI = np.zeros( len(SimNames))
    stdTotDistSI = np.zeros( len(SimNames))
    totDistSG = np.zeros( len(SimNames))
    stdTotDistSG = np.zeros( len(SimNames))
    crossSection = np.array([0., 0.1, 0.3, 1.0])

    for i, iSimName in enumerate(SimNames):
        finalPickleFile = 'Pickles/finalHalos_'+iSimName+'.pkl'

        if os.path.isfile( finalPickleFile ):
            iClusterSample = pkl.load(open(finalPickleFile,'rb'))
        else:
            iClusterSample = ClusterSample(iSimName)

            iClusterSample.extractMergerHalos()

            iClusterSample.CalculateOffsetVectors(nClusters=5)
            
            pkl.dump(iClusterSample, open(finalPickleFile, 'wb'))
        
        
        index = iClusterSample.dist_sg < 250
        
        nClusters = len(index[index])
        print 'rayleigh',rayleigh.fit(iClusterSample.dist_sd[index], loc=0)
        print("selected %f/%f clusters" % \
                  (nClusters,len(iClusterSample.ClusterMass)))

        #plotBulletVector( iClusterSample )

        
        meanSG, stdSG = getMean(iClusterSample.dist_sg[index] )
        meanSI, stdSI = getMean(iClusterSample.dist_si[index])
        meanDI, stdDI = getMean(iClusterSample.dist_di[index])
        meanBeta, stdBeta = getMeanCauchy(iClusterSample.beta[index])
        
        #plotX = np.linspace(-30.,30., 1000)
        #axarr[0].plot(plotX, norm.pdf(plotX, meanSI, stdSI*np.sqrt(nClusters)))
        totDistSI[i] = meanSI
        stdTotDistSI[i] = stdSI
        totDistSG[i] = meanSG
        stdTotDistSG[i] = stdSG
        totBeta[i] = meanBeta
        stdTotBeta[i] = stdBeta
        

        
        print("%s SG mean is %f +/- %0.2f" %  ( iSimName, meanSG, stdSG))
        print("%s SI mean is %f +/- %0.2f" %  ( iSimName, meanSI, stdSI))
        print("%s DI mean is %f +/- %0.2f" %  ( iSimName, meanDI, stdDI))
        print("%s Beta mean is %f +/- %0.2f" %  ( iSimName, meanBeta, stdBeta))
    #totBeta -=  totBeta[0]
    
    #Dump this info in a pickle so i can use it for other scripts
    pkl.dump([crossSection,totBeta,stdTotBeta, totDistSI, totDistSG],\
                 open('Pickles/SimulationBetaToCross.pkl','wb'))
    #plotTrend( crossSection, totDistSI, None, 'Simulations', plotLims=[crossSection[0],3.] )
    #axarr[0].errorbar( crossSection, totDistSI, yerr=stdTotDistSI, fmt='o')
    #axarr[0].set_xlim(1e-2, 2.0)
    #axarr[0].set_xscale('log')
    plt.sca(axarr)
    ax = plt.gca()

    BetaCrossSectionTrend = \
      plotTrend( crossSection, totBeta, stdTotBeta, 'Simulations', \
                   plotLims=[1e-3,3.] )
                   
    ax.errorbar( crossSection, totBeta, yerr=stdBeta, color='black', fmt='o', \
                           label='BAHAMAS-SIDM hydro simulations' )
    crossSectionPlot = 10**np.linspace(-3,np.log10(2),100)
    
    ax.errorbar( crossSectionPlot, 1.-np.exp(-1.*crossSectionPlot/6.5), ls='-', \
                       label='Approximate H14 analytical model', color='k')

                       
    HarveyMeasurement = np.zeros(2)+0.07
    HarveyCrossSection =   betaToCross( HarveyMeasurement, 1., 6.5)
    BAHAMASCrossSection =  betaToCross( HarveyMeasurement, \
                                            *BetaCrossSectionTrend)



    #plotJauzacMeasurement( BetaCrossSectionTrend )
    
    ax.errorbar( [1e-3, BAHAMASCrossSection[0]], \
                           HarveyMeasurement, color='red')
    ax.text( 0.012, HarveyMeasurement[0]*1.1, \
        'Existing HST+Chandra archive of 30 clusters', color='red')
    #ax.errorbar(  HarveyCrossSection, \
    #            [-0.1,HarveyMeasurement[0]], color='red', ls='--', \
    #                        label='Published particle constraints')
    #ax.errorbar(  BAHAMASCrossSection, \
    #            [-0.1,HarveyMeasurement[0]], color='red',ls=':',\
    #                    label='Revised interpretation (cf.W17)')


    plotSensitivity(meanSI,meanSG,BetaCrossSectionTrend)
    

    ax.set_xscale('log')
    ax.set_yscale('log')

    #axarr[0].legend(loc=0)
    ax.legend(loc=0)

    #axarr[0].set_xlabel('Star - Gas Distance / kpc')
    #axarr[0].set_ylabel('Star - Intersection / kpc')
    ax.set_xlabel(r'Dark Matter Cross-Section / cm$^2$/g', fontsize=12)
    ax.set_ylabel(r'Observable $\beta=\delta_{\rm SI}/\delta_{\rm SG}$', fontsize=12)
    ax.set_ylim(2e-3,0.3)
    ax.set_xlim(1e-2,1.75)
    plt.savefig('../plots/HSTconstraints.pdf')
    plt.show()
def getMaxLike( y, bins=30 ):

    #mean = np.median(x)
    #std = np.std(x)/np.sqrt(len(x))

    ypdf, x = np.histogram( y, bins=bins, density=True)
    xc = (x[:-1]+x[1:])/2.
    maxLike = xc[ np.argmax(ypdf) ]
    deltax = xc[1] - xc[0]
    cumsum = np.cumsum(ypdf*deltax)
    std = x[np.argmin(np.abs(cumsum - 0.84))] - maxLike
    
    return maxLike, std/np.sqrt(len(y))


def getMean( x, ax=None ):

    x = np.sort(x)
    mean, std = norm.fit( x )

    

    return np.mean(x), std/np.sqrt(len(x))


def getMeanCauchy( x ):


    mean, std = cauchy.fit( x )
    
    
    return np.mean(x), std/np.sqrt(len(x))


def plot( x, y, label, c):
    nbins = 7
    bins = np.linspace(0,100,nbins+1)

    meany = np.zeros( nbins )
    std = np.zeros( nbins )

    for i in xrange(nbins):

        inbin = ( x > bins[i] ) & (x <= bins[i+1])

        meany[i] = np.mean( y[inbin] )
        std[i] = np.std( y[inbin] ) / np.sqrt(len(y[inbin]))

    bincentres = (bins[:-1] + bins[1:])/2.
    bincentres = bincentres[np.isfinite(std)]
    meany = meany[np.isfinite(std)]
    ax = plt.gca()
    std[:] = np.mean(std)
    plotTrend( bincentres, meany, None, label, color=c )
    ax.errorbar( bincentres, meany, yerr=std, color=c, label=label, fmt='o')
        






def plotSensitivity(si,sg,BetaCrossSectionTrend):
    '''
    Plot the sensitivity in beta
    '''
    error =60.
    sigmaBeta = 0.07
    nClusterBenchMark = 72
    nNewClusters = 72
    totalClusters = nClusterBenchMark + nNewClusters

    ax = plt.gca()
    plotCrossSections = np.linspace(0.001,2.0,2.)
    SensitivityBetaSingleBand = np.zeros(len(plotCrossSections))+sigmaBeta
    SensitivityBetaDualBand = np.zeros(len(plotCrossSections))+sigmaBeta/2.

    SensitivityBetaSingleBand *= np.sqrt(nClusterBenchMark)
    SensitivityBetaDualBand *= np.sqrt(nClusterBenchMark)
    yerr = np.zeros(len(plotCrossSections))+0.01



    BetaLimit = SensitivityBetaSingleBand[0]/np.sqrt(totalClusters)
    CrossSectionLimit =  betaToCross( BetaLimit, *BetaCrossSectionTrend)

    
    #ax.errorbar( [plotCrossSections[0],CrossSectionLimit], \
    #             [BetaLimit,BetaLimit], ls='-', \
    #             color='orange')
    #ax.text( 0.006, BetaLimit, \
    #            '85 Clusters + single band', color='orange')

    #ax.errorbar( [CrossSectionLimit,CrossSectionLimit],\
     #            [-1,BetaLimit], ls='-', \
     #                color='orange')
    
    
    BetaLimit = SensitivityBetaDualBand[0]/np.sqrt(totalClusters)
    CrossSectionLimit =  betaToCross( BetaLimit, *BetaCrossSectionTrend)
      
    ax.errorbar( [plotCrossSections[0],CrossSectionLimit], \
                 [BetaLimit,BetaLimit], ls='-', \
                color='green')
    ax.errorbar( [CrossSectionLimit,CrossSectionLimit], \
                 [1e-4,BetaLimit], ls='-', \
                color='green')
                
    ax.text( 0.012, BetaLimit*1.1, \
        '$\sim10^2$ Clusters Euclid', color='green')

    print 'Cross section limit is ',CrossSectionLimit
    ax.errorbar( [CrossSectionLimit,CrossSectionLimit],\
                 [-1,BetaLimit], ls='-', \
                     color='green')


def plotJauzacMeasurement(BetaCrossSectionTrend):

    dm = 0.5
    dmerr = 0.5
    gas = 1.7
    gaserr = 0.5


    JauzacMeasurement = dm/gas
    JauzacMeasurementError = np.sqrt((dmerr/dm)**2+(gaserr/gas)**2)*JauzacMeasurement
    JauzacCrossSection = (JauzacMeasurement - \
      BetaCrossSectionTrend[0])/BetaCrossSectionTrend[1]
    JauzacCrossSectionError = \
      JauzacMeasurementError/JauzacMeasurement*JauzacCrossSection

    ax = plt.gca()

    ax.plot( [1e-3, JauzacCrossSection], \
                 [JauzacMeasurement,JauzacMeasurement], \
                 '-', color='cyan', label='J18')

    ax.fill_between( [1e-3, 10], \
                    [JauzacMeasurement-JauzacMeasurementError,\
                         JauzacMeasurement-JauzacMeasurementError], \
                    [JauzacMeasurement+JauzacMeasurementError,\
                         JauzacMeasurement+JauzacMeasurementError], \
                         color='cyan', alpha=0.5)
                         
    ax.fill_between( [JauzacCrossSection-JauzacCrossSectionError, \
                          JauzacCrossSection+JauzacCrossSectionError], \
                    [JauzacMeasurement,JauzacMeasurement], \
                    [-1,-1], \
                         color='cyan', alpha=0.5)
    print("JAUZAC CONSTRAINTS ARE %0.2f +/- %0.2f" \
              %(JauzacCrossSection,JauzacCrossSectionError))

    x=np.linspace(-10,30,1000)
    gaussCDF = np.cumsum(norm.pdf(x, JauzacCrossSection, JauzacCrossSectionError))
    gaussCDF/=np.max(gaussCDF)
    JauzacCrossSectionLimit = x[np.argmin(np.abs(gaussCDF-0.68))]
    print("JAUZAC LIMITS IS <  %0.2f " \
              %(JauzacCrossSectionLimit))

    
    
def betaToCross( beta, Aparam, sigmaStar):
    return np.log(1.-beta/Aparam)*-sigmaStar


def plotBulletVectors( iClusterSample ):
     plt.sca( axarr[0] )
     plot( iClusterSample.dist_sg[index], iClusterSample.dist_si[index], \
                  iSimName, colors[i])

     axarr.scatter( iClusterSample.dist_sg[index], iClusterSample.dist_si[index], \
                          color=colors[i])

     axarr.hist(  iClusterSample.dist_sd[index], bins=30, alpha=1-i*0.3, \
                           density=True, color=colors[i] )
     axarr[1].hist(  iClusterSample.dist_sg[index], bins=30, alpha=1-i*0.3,\
                            density=True)
     axarr[2].hist(  iClusterSample.beta[index], bins=30, alpha=1-i*0.3,\
                            density=True)
if __name__ == '__main__':
    HSTconstraints()
