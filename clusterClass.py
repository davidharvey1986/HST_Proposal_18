import component_extractor as ce
import pyfits as fits
import numpy as np
import os as os
import matchMassComponents as mmc
'''
clusterClass just gets each cluster and source extracts all the halos
from each mass component.

'''


class clusterClass:
    
    def __init__( self, clusterInt, simulation, \
                      redshift, GasConcentration ):
        redshiftStr = "%0.3f" % redshift

        self.initAllPositions()
        self.GasConcentration = \
          GasConcentration

        self.dataDir = \
          '/Users/DavidHarvey/Documents/Work/Mergers/sims/BAHAMAS'
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

        self.getClusterMass()
        if self.GasConcentration < 0.2:
            self.getClusterMembers()
            self.getGasPosition()
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
        
    def getClusterMass( self ):
        dtypes = [('id', object), ('fof', int), \
                      ('Mass', float), ('R200', float)]

        AllClusters = \
          np.loadtxt( self.clusterCatalog, dtype=dtypes)

        self.mass = \
          AllClusters['Mass'][ AllClusters['id'] == self.clusterInt ][0]
          

    def getClusterMembers( self ):

        dtypes =  [('id', object), ('x', float), \
                      ('y', float), ('z', float), \
                       ('mass30', float), ('mass100', float)]


        self.clusterMembers = \
          np.loadtxt( self.GalaxyCat, dtype=dtypes)

        
    def getGasPosition( self, smoothing=5.):
        #Not the gas smoothing scale is x5kpc not x1kpc
        allGasSources = \
          ce.component_extractor( self.GasImage, smoothing=smoothing, \
                                redshift=self.redshift,\
                                    filtername='gauss_5.0_9x9.conv')
        #The xray maps have a different resolution to that
        #of the other maps (5kpc, as supposed  from 1kpc)
        try:
            self.xGasPositions = (allGasSources['X_IMAGE'] - 1000.)*5. + 1000
            self.yGasPositions = (allGasSources['Y_IMAGE'] - 1000.)*5. + 1000
        except:
            self.xGasPositions = np.array([1000])
            self.yGasPositions = np.array([1000])
        
    def getStellarPosition( self ):

        allStellarSources = \
          ce.component_extractor( self.StellarImage, smoothing=0., \
                                redshift=self.redshift,\
                                filtername='gauss_5.0_9x9.conv')
                                
        self.xStellarPositions = allStellarSources['X_IMAGE']
        self.yStellarPositions = allStellarSources['Y_IMAGE']


    def getDarkMatterPositions( self ):

        stellarMass = fits.open( self.StellarImage )[0].data
        totalMass = fits.open( self.TotalImage )[0].data
        dmMass = totalMass - stellarMass
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
        
        
    def getTotalPosition( self ):

        allTotalSources = \
          ce.component_extractor( self.TotalImage, smoothing=0., \
                                redshift=self.redshift,\
                                filtername='gauss_5.0_9x9.conv')
                                
        self.xTotalPositions = allTotalSources['X_IMAGE']
        self.yTotalPositions = allTotalSources['Y_IMAGE']
    

    def checkCluster( self ):
         #Check that there are found gas posisionts
        if type(self.xGasPositions) == float:
            raise ValueError('Cannot match since no gas positions available')

        #CHeck if those positions are within the field of view of the
        #hires shit
        GasRadialDistance  = \
          np.sqrt( ( self.xGasPositions - 1000)**2 + \
                   ( self.yGasPositions - 1000)**2)

        if len( GasRadialDistance[ GasRadialDistance < 1000 ]) < 2:
            print GasRadialDistance, self.clusterInt
            raise ValueError('No gas positions in DM FOV available')           
    
    def combineMassComponents( self ):
        #Combine the three mass components
        self.checkCluster()
        
        GasDarkMatter = \
          mmc.matchMassComponents( self.xGasPositions, \
                                   self.yGasPositions, \
                                   self.xDarkMatterPositions, \
                                   self.yDarkMatterPositions, \
                                   searchRadKpc = 100.)

        GasStellar = \
          mmc.matchMassComponents( self.xGasPositions, \
                                   self.yGasPositions, \
                                   self.xDarkMatterPositions, \
                                   self.yDarkMatterPositions, \
                                   searchRadKpc = 100.)
