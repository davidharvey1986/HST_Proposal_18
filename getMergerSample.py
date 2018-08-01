
from clusterClass import *
import getSimData as gSD
import os as os
import glob as glob
import pickle as pkl
import xray_sim_concentration as xsc
from MockedClusterClass import *

class ClusterSample:

    def __init__( self, SimName, **kwargs ):
        self.SimName = SimName
        self.ClusterSample = []
        if 'Mocked' in SimName:
            self.getMockedSample( **kwargs )
        
        else:
            self.Sims = gSD.getSimList( SimName )
            self.dataDir = \
            '/Users/DavidHarvey/Documents/Work/Mergers/sims/BAHAMAS'
            self.pickleFile = 'Pickles/'+SimName+'_MassPositions.pkl'

            if os.path.isfile( self.pickleFile):
                self.ClusterSample = \
                pkl.load( open( self.pickleFile, 'rb'))
            else:
                self.getClusterSample()

    def getClusterSample( self ):
        
        for iSim in self.Sims:
          
                self.SimRoot = self.dataDir+'/'+iSim
                self.getRedshiftList( iSim )

                for iRedshift in self.redshiftList:
                    simPickleFile = \
                      "Pickles/"+iSim+"_"+str(iRedshift)+".pkl"

                    if os.path.isfile( simPickleFile ):
                        simClusterClasses = \
                          pkl.load( open(simPickleFile, 'rb'))
                        [ self.ClusterSample.append( i ) for i in \
                        simClusterClasses ]
                        
                    else:
                        self.getClusterInts( iRedshift )
                    
                        self.getGasConcentration( iSim, iRedshift)
                
                        for iCluster, iClusterInt in \
                          enumerate(self.clusterInts):
                            print 'ClusterInt',iClusterInt
                            iClusterClass = \
                            clusterClass( iClusterInt, iSim, iRedshift,\
                                    self.GasConcentration[iCluster])
                        
                            self.ClusterSample.append( iClusterClass )
                        
                        pkl.dump( self.ClusterSample, open(simPickleFile,'wb'))

        pkl.dump( self.ClusterSample, open(self.pickleFile, 'wb') )
                    
                                                 
            

    def getGasConcentration( self, sim, redshift ):
        redshiftStr = "z_%0.3f" % redshift
        self.GasConcentration = \
          xsc.get_all_xray_concentrations( self.clusterInts, \
                                        sim=sim,\
                                        clusterRedshiftStr=redshiftStr )
        

        
    def getRedshiftList( self, iSim ):

        redshiftListComplete = glob.glob(self.SimRoot+'/z_*')
        print redshiftListComplete
        
        self.redshiftList = [ np.float(i.split('/')[-1].split('_')[1])\
          for i in redshiftListComplete]


        
    def getClusterInts( self, redshift ):
        redshiftStr = "%0.3f" % redshift
        self.GalaxyCatalogs = \
          glob.glob(self.SimRoot+'/z_'+redshiftStr+'/GALAXY_CATALOGS/cluster*')

        self.clusterInts = [ i.split('/')[-1].split('_')[1] \
                                 for i in self.GalaxyCatalogs ]

        
    def extractMergerHalos( self ):

        mergerHalos = []
        pickleFile = "Pickles/"+self.SimName+"_mergerHalos.pkl"
        if os.path.isfile( pickleFile ):
            self.mergerHalos = pkl.load(open(pickleFile, 'rb'))
        else:
            
            for iCluster in self.ClusterSample:

                flag = iCluster.combineMassComponents()
                if flag > 0:
                    iCluster.calculateBulletVectors()
                    mergerHalos.append( iCluster )

                
                    
            self.mergerHalos = mergerHalos
            pkl.dump( mergerHalos, open(pickleFile, 'wb'))

    def CalculateOffsetVectors( self, nClusters=3 ):

        self.dist_si = np.array([])
        self.dist_sg = np.array([])
        self.dist_di = np.array([])
        self.beta = np.array([])
        self.betaPerp = np.array([])
        self.ClusterMass = np.array([])
        self.dist_sd = np.array([])
        for iMerger in self.mergerHalos:
                
            if (len(iMerger.dist_si)<=1) | (len(iMerger.dist_si) > nClusters):
                continue
            self.dist_si = np.append( self.dist_si, iMerger.dist_si[0:nClusters])
            self.dist_sg = np.append( self.dist_sg, iMerger.dist_sg[0:nClusters])
            self.dist_di = np.append( self.dist_di, iMerger.dist_di[0:nClusters])
            sd = np.sqrt(np.sum(iMerger.vector_sd**2, axis=0))
            self.dist_sd = np.append( self.dist_sd, sd[0:nClusters])

            self.beta = np.append( self.beta, iMerger.beta[0:nClusters])
            self.betaPerp = np.append( self.betaPerp, iMerger.betaPerp[0:nClusters])
            self.ClusterMass = \
              np.append( self.ClusterMass, \
                             np.zeros(len(iMerger.beta[0:nClusters]))+iMerger.mass)
            

            


    def getMockedSample( self, **kwargs ):

        for iCluster in xrange(kwargs['nClusters']):
            iClusterSample = \
            MockedCluster( kwargs['DarkMatterStellarOffset'],  \
                    kwargs['GasStellarOffset'],
                    kwargs['DarkMatterError'], \
                    kwargs['GasError'],\
                    kwargs['nSpuriousHalos'])
            self.ClusterSample.append( iClusterSample )
        

    
            
