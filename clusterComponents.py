'''
This will loop through the all the sims and get all the stellar and dark matter
positions
'''

import component_extractor as ce
import glob
import numpy as np
import pyfits as fits
import os
import getSimData as gsd
import sys

def main( simName ):
    
    
    simulationList = gsd.getSimList( simName )
    for iSimulation in simulationList:
        getPositionsForGivenSimulation(iSimulation)
        
def getPositionsForGivenSimulation( simulation ):
    '''
                      get the stellar and dark amtter positions for a given 
    simulation

    '''
    dataDir="/data1/BAHAMAS"
    redshiftList = getRedshiftList( dataDir+'/'+simulation)
    for iRedshift in redshiftList:
        clusterInts = getClusterInts( dataDir+'/'+simulation, \
                                     iRedshift )
        redshiftStr = "%0.3f" % iRedshift


        gasConcentration = \
            get_all_xray_concentrations( clusterInts, sim=simulation, \
                                   clusterRedshiftStr='z_'+redshiftStr)
        
        listOfClusterComponents = []
        
        for iGas, iClusterInt in enumerate(clusterInts):
            redshiftStr = "%0.3f" % iRedshift
            
            pklFile = "Pickles/StellarAndDM_%s_z_%s_%s.pkl" %(simulation,iRedshiftStr,iClusterInt)
            if os.path.isfile( pklFile ):
                print("File already done")
                continue
            
            iClusterComponent = \
                clusterComponents( iClusterInt, simulation, \
                                   iRedshift, iGas)
            
            pkl.dump(iClusterComponent, open(pklFile,'wb'))
            
            
def getClusterInts( sim, redshift):
    redshiftStr = "%0.3f" % redshift

    GalaxyCatalogs = \
        glob.glob(sim+'/z_'+redshiftStr+'/GALAXY_CATALOGS/cluster*')

    clusterInts = [ i.split('/')[-1].split('_')[1] \
                         for i in GalaxyCatalogs ]
    return clusterInts
        
def getRedshiftList( iSim ):

    redshiftListComplete = glob.glob(iSim+'/z_*')
    
    redshiftList = [ np.float(i.split('/')[-1].split('_')[1])\
                          for i in redshiftListComplete]
   
    return redshiftList

class clusterComponents:
    
    def __init__( self, clusterInt, simulation, \
                      redshift, GasConcentration ):
        redshiftStr = "%0.3f" % redshift
        self.initAllPositions()
        self.GasConcentration = \
          GasConcentration

        self.dataDir = \
          '/data1/BAHAMAS/'
        self.clusterInt = clusterInt
        self.simulation = simulation
        self.redshift = redshift
        
        self.SimRoot = self.dataDir+'/'+simulation+'/z_'+str(redshiftStr)
        
        self.GasImage = \
          self.SimRoot+'/MAPS/cluster_'+str(clusterInt)+'_xray.fits'
          
        self.StellarImage = \
          self.SimRoot+'/HIRES_MAPS/cluster_'+str(clusterInt)+'_stellar_sph.fits'

        self.TotalImage = \
          self.SimRoot+'/HIRES_MAPS/cluster_'+str(clusterInt)+'_total_sph.fits'
          
        self.GalaxyCat = \
          self.SimRoot+'/GALAXY_CATALOGS/cluster_'+str(clusterInt)+\
          '_galaxy_catalog.dat'

        self.clusterCatalog = self.SimRoot+'/catalog.dat'


        if self.GasConcentration < 0.2:
            self.getStellarPosition()
            self.getDarkMatterPositions()
    def initAllPositions( self ):
        self.xGasPositions = 1000.
        self.yGasPositions = 1000.
        self.xStellarPositions = 1000.
	self.yStellarPositions = 1000.
        self.xTotalPositions = 1000.
        self.yTotalPositions = 1000.
        self.xDarkMatterPositions = 1000.
	self.yDarkMatterPositions = 1000.

       
    
    def getStellarPosition( self ):

        allStellarSources = \
          ce.component_extractor( self.StellarImage, smoothing=0., \
                                redshift=self.redshift,\
                                filtername='gauss_5.0_9x9.conv')
        
        if 'X_IMAGE' in allStellarSources.dtype.names:
            self.xStellarPositions = allStellarSources['X_IMAGE']
            self.yStellarPositions = allStellarSources['Y_IMAGE']


    def getDarkMatterPositions( self ):

        stellarMass = fits.open( self.StellarImage )[0].data
        totalMass = fits.open( self.TotalImage )[0].data
        gasMass = fits.open( self.GasImage )[0].data
        dmMass = totalMass - stellarMass - gasMass
        
        randomStr = str(np.random.random_integers(0,10000))
        DMfile = 'DMmass_'+str(randomStr)+'.fits'
        fits.writeto( DMfile, dmMass )
        
        allDMSources = \
          ce.component_extractor( DMfile, smoothing=0., \
                                redshift=self.redshift,\
                                filtername='gauss_5.0_9x9.conv')
                                
        self.xDarkMatterPositions = allDMSources['X_IMAGE']
        self.yDarkMatterPositions = allDMSources['Y_IMAGE']
        os.system('rm -fr '+DMfile)
    

    

if __name__ == '__main__':

    main(sys.argv[1])
