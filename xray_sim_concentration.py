import pickle as pkl
import numpy as np
def get_all_xray_concentrations( clusterInts, \
                                     sim='AGN_TUNED_nu0_L400N1024_WMAP9',\
                                     clusterRedshiftStr='z_0.250' ):
    '''
    clusterInts : the integers of each cluster name
    loop through all xray maps and get concentrations
    '''

    dataDir = '/Users/DavidHarvey/Documents/Work/Mergers/code/MergersBAHAMAS/Pickles/GalaxyCataloguesPickles'

    pklFile = dataDir+'/'+sim+'_'+clusterRedshiftStr+'.pkl'

    refClusterInts, Xray, offset = \
      pkl.load(open(pklFile,'rb'))
    xray_concentration = []
    for i in clusterInts:
        iXray = np.array(Xray)[ i == refClusterInts.astype(str) ]
        if len(iXray) == 0:
            xray_concentration.append()
        else:
            xray_concentration.append(iXray[0])

    
   
    return np.array(xray_concentration)
    
