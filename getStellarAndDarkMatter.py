from clusterClass import *
'''
This is the class used to run sextractor over all
the clusters on Ketel to speed it up
'''


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
