from getMergerSample import *
import getSimData as GSD
from matplotlib import pyplot as plt

from scipy.stats import norm
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import cauchy
from scipy.stats import rayleigh

def HSTconstraints():

    '''
    plot and devise the constraints for the HST proposal
    '''
    SimNames = ['CDM', 'SIDM0.1','SIDM'] #GSD.getSimNameList()
    fig, axarr = plt.subplots( 1, figsize=(8,20) )

    bins=40
    colors = ['k','r','b']

    totBeta = np.zeros( len(SimNames))
    stdTotBeta = np.zeros( len(SimNames))
    crossSection = np.array([0., 0.1, 1.0])
    
    for i, iSimName in enumerate(SimNames):

        iClusterSample = ClusterSample(iSimName)

        iClusterSample.extractMergerHalos()

        iClusterSample.CalculateOffsetVectors(nClusters=1)
        

    
        index = np.abs(iClusterSample.dist_si) > 0.
        
        
        nClusters = len(index[index])
        print 'rayleigh',rayleigh.fit(iClusterSample.dist_sd[index], loc=0)
        print("selected %f/%f clusters" % \
                  (nClusters,len(iClusterSample.ClusterMass)))
        
        #plot( iClusterSample.dist_sg[index], iClusterSample.beta[index], \
        #          iSimName, colors[i])
        #axarr.scatter( iClusterSample.dist_sg[index], iClusterSample.dist_si[index], \
        #                   color=colors[i])

        axarr.hist(  iClusterSample.dist_sd[index], bins=30, alpha=1-i*0.3, \
                            density=True )
        #axarr[1].hist(  iClusterSample.dist_sg[index], bins=30, alpha=1-i*0.3,\
        #                    density=True)
        #axarr[2].hist(  iClusterSample.beta[index], bins=30, alpha=1-i*0.3,\
        #                    density=True)

        
        meanSG, stdSG = getMean(iClusterSample.dist_sg[index] )
        meanSI, stdSI = getMean(iClusterSample.dist_si[index])
        meanDI, stdDI = getMean(iClusterSample.dist_di[index])
        meanBeta, stdBeta = getMeanCauchy(iClusterSample.beta[index])
        
        #plotX = np.linspace(-30.,30., 1000)
        #axarr[0].plot(plotX, norm.pdf(plotX, meanSI, stdSI*np.sqrt(nClusters)))
        
        totBeta[i] = meanBeta
        stdTotBeta[i] = stdBeta
        
        
        
        print("%s SG mean is %f +/- %0.2f" %  ( iSimName, meanSG, stdSG))
        print("%s SI mean is %f +/- %0.2f" %  ( iSimName, meanSI, stdSI))
        print("%s DI mean is %f +/- %0.2f" %  ( iSimName, meanDI, stdDI))
        print("%s Beta mean is %f +/- %0.2f" %  ( iSimName, meanBeta, stdBeta))
        
    #axarr.legend()
    #axarr.errorbar( crossSection, totBeta, yerr=stdBeta, fmt='o' )
    #axarr.plot( crossSection, 1.-np.exp(-1.*crossSection/6.5), '-')
    
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


    mean, std = norm.fit( x )

    return np.mean(x), np.std(x)/np.sqrt(len(x))


def getMeanCauchy( x ):


    mean, std = cauchy.fit( x )
    
    return np.mean(x), np.std(x)/np.sqrt(len(x))


def plot( x, y, label, c):
    nbins = 10
    bins = np.linspace(0,100,nbins+1)

    meany = np.zeros( nbins )
    std = np.zeros( nbins )

    for i in xrange(nbins):

        inbin = ( x > bins[i] ) & (x <= bins[i+1])

        meany[i] = np.mean( y[inbin] )
        std[i] = np.std( y[inbin] ) / np.sqrt(len(y[inbin]))
    ax = plt.gca()
    bincentres = (bins[:-1] + bins[1:])/2.
    ax.errorbar( bincentres, meany, yerr=std, color=c, label=label, ls='-')
        



