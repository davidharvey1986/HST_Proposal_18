from getMergerSample import *
import getSimData as GSD
from matplotlib import pyplot as plt

from scipy.stats import rayleigh
from scipy.stats import kurtosis

def HSTconstraints():

    '''
    plot and devise the constraints for the HST proposal
    '''
    SimNames = ['CDM']#, 'SIDM0.1','SIDM'] #GSD.getSimNameList()
    fig, axarr = plt.subplots( 1, figsize=(8,8) )
    bins=40
    colors = ['k','r','b']
    for i, iSimName in enumerate(SimNames):

        iClusterSample = ClusterSample(iSimName)

        iClusterSample.extractMergerHalos()

        iClusterSample.CalculateOffsetVectors(nClusters=3)
        

    
        index = np.abs(iClusterSample.dist_si) < 20.
        
        
  
        #p#rint 'rayleigh',rayleigh.fit(iClusterSample.dist_sd[index], loc=0)
        print("selected %f/%f clusters" % \
                  (len(index[index]),len(iClusterSample.ClusterMass)))
        
        plot( iClusterSample.dist_sg[index], iClusterSample.dist_si[index], iSimName, colors[i])
        axarr.scatter( iClusterSample.dist_sg[index], iClusterSample.dist_di[index], color=colors[i])

        #axarr[0].hist(  iClusterSample.dist_si[index], bins=30 )
        #axarr[1].hist(  iClusterSample.dist_sg[index], bins=30 )
        #axarr[2].hist(  iClusterSample.beta[index], bins=np.linspace(-1,1,41) )

        
        meanSG, stdSG = getMean(iClusterSample.dist_sg[index] )
        meanSI, stdSI = getMean(iClusterSample.dist_si[index])
        meanDI, stdDI = getMean(iClusterSample.dist_di[index])
        meanBeta, stdBeta = getMean(iClusterSample.beta[index])
        
        

        print("%s SG mean is %f +/- %0.2f" %  ( iSimName, meanSG, stdSG))
        print("%s SI mean is %f +/- %0.2f" %  ( iSimName, meanSI, stdSI))
        print("%s DI mean is %f +/- %0.2f" %  ( iSimName, meanDI, stdDI))
        print("%s Beta mean is %f +/- %0.2f" %  ( iSimName, meanBeta, stdBeta))
        
    #axarr.legend()
        
def getMaxLike( y, bins ):

    #mean = np.median(x)
    #std = np.std(x)/np.sqrt(len(x))

    ypdf, x = np.histogram( y, bins=bins, density=True)
    xc = (x[:-1]+x[1:])/2.
    maxLike = xc[ np.argmax(ypdf) ]
    deltax = xc[1] - xc[0]
    cumsum = np.cumsum(ypdf*deltax)
    std = x[np.argmin(np.abs(cumsum - 0.84))] - maxLike
    
    return maxLike, std/np.sqrt(len(y))


def getMean( x ):

    return np.mean(x), np.std(x)/np.sqrt(len(x))


def plot( x, y, label, c):
    nbins = 6
    bins = np.linspace(0,100,nbins+1)

    meany = np.zeros( nbins )
    for i in xrange(nbins):

        inbin = ( x > bins[i] ) & (x <= bins[i+1])

        meany[i] = np.mean( y[inbin] )


    ax = plt.gca()
    bincentres = (bins[:-1] + bins[1:])/2.
    ax.plot( bincentres, meany, c, label=label, ls='-')
        
