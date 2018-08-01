from getMergerSample import *
from matplotlib import pyplot as plt
import itertools
import ipdb as pdb
def checkMergerMatching( SimName ):

    '''
    CDM seems to have some massive bias
    Need to check out what is going on
    
    First check the matching

    '''
    marker = itertools.cycle(('D', '+','<','*')) 

    Clusters = ClusterSample(SimName)

    Clusters.extractMergerHalos()
    fig = plt.gcf()
    for iCluster in Clusters.mergerHalos:

        plt.clf()
        plt.xlim(0,2000)
        plt.ylim(0,2000)
        plt.plot( iCluster.xGasPositions, \
                  iCluster.yGasPositions, 'r.')
        plt.plot( iCluster.xDarkMatterPositions, \
                  iCluster.yDarkMatterPositions, 'b.')
        plt.plot( iCluster.xStellarPositions, \
                  iCluster.yStellarPositions, 'g.')

        for i, iMerger in  enumerate(iCluster.mergerHalos):
            iMarker = marker.next()
            plt.plot( iMerger['xGas'], iMerger['yGas'], \
                          marker=iMarker, color='red' )
            plt.plot( iMerger['xDarkMatter'], iMerger['yDarkMatter'], 
                          marker=iMarker, color='blue' )
            plt.plot( iMerger['xStellar'], iMerger['yStellar'], \
                          marker=iMarker, color='green' )
        print("There are %i mergers in cluster %s" % \
                  (len(iCluster.mergerHalos), iCluster.clusterInt))
        plt.draw()
        fig.canvas.flush_events()
        raw_input("Press Enter to continue...")
        
 
def checkMockedClusters():
    fig, axarr = plt.subplots( 5, figsize=(4,20) )

    SimName = 'MockedCDM'
    mergerSample = \
      ClusterSample('MockedCDM', nClusters=100, \
                        DarkMatterStellarOffset=0., \
                        GasStellarOffset=50., \
                        DarkMatterError=10., \
                        GasError=10.,\
                        nSpuriousHalos=20)

    mergerSample.extractMergerHalos()
    mergerSample.CalculateOffsetVectors()
    index = (mergerSample.ClusterMass >  14.)
    
    axarr[0].hist( mergerSample.dist_si[index], \
                           bins=50, label=SimName )
    axarr[1].hist( mergerSample.dist_sg[index], bins=50 )
    axarr[2].hist( mergerSample.dist_di[index], bins=50 )
    axarr[3].hist( mergerSample.betaPerp[index], bins=100 )
    axarr[4].hist( mergerSample.beta[index], bins=100 )
    
    meanSG, stdSG = getMean(mergerSample.dist_sg[index] )
    meanSI, stdSI = getMean(mergerSample.dist_si[index])
    meanDI, stdDI = getMean(mergerSample.dist_di[index])
    meanBeta, stdBeta = getMean(mergerSample.beta[index])

        
    print("%s SG mean is %f +/- %0.2f" %  ( SimName, meanSG, stdSG))
    print("%s SI mean is %f +/- %0.2f" %  ( SimName, meanSI, stdSI))
    print("%s DI mean is %f +/- %0.2f" %  ( SimName, meanDI, stdDI))
    print("%s Beta mean is %f +/- %0.2f" %  ( SimName, meanBeta, stdBeta))
        
    axarr[0].legend()
def getMean( x ):

    mean = np.mean(x)
    std = np.std(x)/np.sqrt(len(x))
    return mean, std
