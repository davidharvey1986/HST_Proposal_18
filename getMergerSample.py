
from clusterClass import *
import getSimData as gSD
import os as os
import glob as glob
import pickle as pkl
import xray_sim_concentration as xsc


class ClusterSample:

    def __init__( self, SimName ):
        self.SimName = SimName
        self.ClusterSample = []
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
                
                self.getClusterInts( iRedshift )
                
                self.getGasConcentration( iSim, iRedshift)
                
                for iCluster, iClusterInt in enumerate(self.clusterInts):
                    iClusterClass = \
                      clusterClass( iClusterInt, iSim, iRedshift,\
                                    self.GasConcentration[iCluster])
                      
                    self.ClusterSample.append( iClusterClass )


        pkl.dump( self.ClusterSample, open(self.pickleFile, 'wb') )
                    
                                                 
            

    def getGasConcentration( self, sim, redshift ):
        redshiftStr = "z_%0.3f" % redshift
        self.GasConcentration = \
          xsc.get_all_xray_concentrations( self.clusterInts, \
                                        sim=sim,\
                                        clusterRedshiftStr=redshiftStr )
        

        
    def getRedshiftList( self, iSim ):

        redshiftListComplete = glob.glob(self.SimRoot+'/z_*')
        self.redshiftList = [ np.float(i.split('/')[-1].split('_')[1])\
          for i in redshiftListComplete]


        
    def getClusterInts( self, redshift ):
        redshiftStr = "%0.3f" % redshift
        self.GalaxyCatalogs = \
          glob.glob(self.SimRoot+'/z_'+redshiftStr+'/GALAXY_CATALOGS/cluster*')[0:20]

        self.clusterInts = [ i.split('/')[-1].split('_')[1] \
                                 for i in self.GalaxyCatalogs ]

        


