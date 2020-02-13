

from HSTconstraints import *

from sklearn.linear_model import LinearRegression

def main():

    nMassBins = 10
    nCross = 10
    SimNames = ['CDM','SIDM0.1','SIDM0.3','SIDM1'] #GSD.getSimNameList()


    allSIforRegression = np.array([])
    regressionFeatureCross  = np.array([])
    regressionFeaturesMass = np.array([])
    for i, iSimName in enumerate(SimNames):
        print iSimName

        finalPickleFile = 'Pickles/finalHalos_'+iSimName+'.pkl'

        if os.path.isfile( finalPickleFile ):
            iClusterSample = pkl.load(open(finalPickleFile,'rb'))
        else:
            iClusterSample = ClusterSample(iSimName)

            iClusterSample.extractMergerHalos()

            iClusterSample.CalculateOffsetVectors(nClusters=5)
            
            pkl.dump(iClusterSample, open(finalPickleFile, 'wb'))
        

        
        massFeature, starIntersectionPoint = \
          getSIfunctionMass( iClusterSample,\
                                 nMassBins=nMassBins )

        if iSimName == 'CDM':
            crossFeature = np.zeros( len(massFeature))
        else:
            crossFeature = np.zeros( len(massFeature))+\
            np.float( iSimName.split('M')[1])

        regressionFeatureCross = np.append(regressionFeatureCross, \
                                            crossFeature)
        regressionFeaturesMass = np.append(regressionFeaturesMass,\
                                            massFeature)
        
        
        allSIforRegression = np.append( allSIforRegression,\
                                    starIntersectionPoint[0,:] )


    allFeatures = np.vstack((regressionFeatureCross,regressionFeaturesMass)).T
    reg = LinearRegression().fit(allFeatures,allSIforRegression)

    dumpForERC = '/Users/DavidHarvey/Dropbox/ERC/code/mergerModel.pkl'
    pkl.dump(reg,open(dumpForERC,'wb'))

    crossTest = 10**np.linspace(-3,0.,nCross)
    crossGrid, massGrid = np.meshgrid( crossTest, massFeature)
    massVector = massGrid.reshape(nCross*nMassBins)
    crossVector = crossGrid.reshape(nCross*nMassBins)
    predictFeatures = np.vstack(( crossVector, massVector)).T

    offsets = reg.predict( predictFeatures)
    
    plots = [ plt.plot( massFeature, offsets[i==crossVector], label=str(i)) \
                  for i in crossTest]
    plt.legend()
    plt.show()
    
def getSIfunctionMass( iClusterSample, nMassBins=5):

    starIntersectionPoint = np.zeros((2,nMassBins), float)
    
    massBins = np.linspace(14.,15.,nMassBins+1)
    for iMassBinIndex in xrange(nMassBins):
        index = (np.abs(iClusterSample.ClusterMass) >  \
                     massBins[iMassBinIndex]) &    \
                     (np.abs(iClusterSample.ClusterMass) < \
                          massBins[iMassBinIndex+1])     
        
        nClusters = len(index[index])

       
        if len(iClusterSample.beta[index]) == 0:
            continue

        distSD = \
          np.sqrt(removeSmoothing(np.abs(iClusterSample.dist_si[index]))**2+\
                      removeSmoothing(np.abs(iClusterSample.dist_di[index]))**2)

        
        meanSD, stdSD = getMean( distSD )

        starIntersectionPoint[0, iMassBinIndex] = np.median(distSD)
          
        starIntersectionPoint[1, iMassBinIndex] = stdSD
          
    massBinCentres = (massBins[1:] + massBins[:-1])/2.
    plt.plot(massBinCentres, starIntersectionPoint[0,:], \
                 label=iClusterSample.SimName)
    return massBinCentres, starIntersectionPoint

def removeSmoothing( dist ):
    '''remove effect of smoothing
    '''
    gamma = 10**(-0.02)
    epsilon = 4.
    alpha = 0.41
    remove =  (np.abs(dist**(gamma) - (alpha*epsilon)))**(1./gamma)

    return remove

def getColor( SimName):
    if SimName == 'SIDM1':
	return 'red'
    if SimName == 'SIDM0.1':
	return 'blue'
    if SimName == 'CDM':
        return 'black'
    if SimName == 'SIDM0.3':
        return 'green'
    if SimName == 'obs':
        return 'cyan'
    if SimName == 'CDMhires':
        return 'black'
    if SimName == 'SIDM1hires':
        return 'red'

if __name__ == '__main__':

    main()
