import component_extractor as ce
import pyfits as fits
import numpy as np
import os as os
import matchMassComponents as mmc
import ipdb as pdb

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

        
    def getGasPosition( self, smoothing=5., rebin_pix=5., \
                            scales='4, 8, 16, 32'):
        #Not the gas smoothing scale is x5kpc not x1kpc
        #allGasSources = \
        #  ce.component_extractor( self.GasImage, smoothing=smoothing, \
        #                        redshift=self.redshift,\
        #                            filtername='gauss_5.0_9x9.conv')

        allGasSources = \
          ce.wavdetectExtractor( self.GasImage, \
                                     scales=scales, \
                                     rebin_pix=rebin_pix)

        #The xray maps have a different resolution to that
        #of the other maps (5kpc, as supposed  from 1kpc)
        dim = 1000./rebin_pix
        scale = 5.*rebin_pix
  
        self.xGasPositions = \
          (allGasSources['X'] - dim)*scale + 1000.
        self.yGasPositions = \
          (allGasSources['Y'] - dim)*scale + 1000.
  
        
    def getStellarPosition( self ):

        allStellarSources = \
          ce.component_extractor( self.StellarImage, smoothing=100., \
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
          ce.component_extractor( DMfile, smoothing=100., \
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
            #print('Cannot match since no gas positions available')
            return -2

        #CHeck if those positions are within the field of view of the
        #hires shit
        GasRadialDistance  = \
          np.sqrt( ( self.xGasPositions - 1000.)**2 + \
                   ( self.yGasPositions - 1000.)**2)

        if len( GasRadialDistance[ GasRadialDistance < 1000 ]) < 2:
            print('No gas positions in DM FOV available %s' % self.clusterInt)
            return -1
        else:
            return 1
    def combineMassComponents( self ):
        #Combine the three mass components
        flag = self.checkCluster()
        if  flag < 0:
            return flag
        
        self.GasDarkMatter = \
          mmc.matchMassComponents( self.xGasPositions, \
                                   self.yGasPositions, \
                                   self.xDarkMatterPositions, \
                                   self.yDarkMatterPositions, \
                                   searchRadKpc = 200.)

        self.GasStellar = \
          mmc.matchMassComponents( self.xGasPositions, \
                                   self.yGasPositions, \
                                   self.xStellarPositions, \
                                   self.yStellarPositions, \
                                   searchRadKpc = 200.)



        
        
        if len(self.GasDarkMatter) == 0 & len(self.GasStellar) == 0:
            return -3                       
        if len(self.GasDarkMatter) != len(self.GasStellar):
            self.RemoveUnmatchedHalos()
        if len(self.GasDarkMatter) != len(self.GasStellar):
            raise ValueError("the remove unmatched didnt work and ive fucked up")
        distSDFlag = self.CheckDarkMatterStellarDistance()
        if distSDFlag < 0:
            return distSDFlag
        

        xGasColumn = \
          fits.Column( name='xGas', \
                           array=self.GasDarkMatter['X_1'], \
                           format='K')
                            
        yGasColumn = \
          fits.Column(  name='yGas', \
                            array=self.GasDarkMatter['Y_1'], \
                            format='K')
                            
        xDarkMatterColumn = \
               fits.Column(name='xDarkMatter', \
                            array=self.GasDarkMatter['X_2'], \
                            format='K')
        yDarkMatterColumn = \
               fits.Column(name='yDarkMatter', \
                           array=self.GasDarkMatter['Y_2'], \
                           format='K')

        xStellarColumn = \
          fits.Column( name='xStellar', \
                        array=self.GasStellar['X_2'], \
                        format='K')
                        
        yStellarColumn = \
          fits.Column(name='yStellar', \
                       array=self.GasStellar['Y_2'], \
                       format='K') 
        
        totalColumns = [xGasColumn, yGasColumn, \
                    xDarkMatterColumn, yDarkMatterColumn, \
                    xStellarColumn, yStellarColumn ]

        
        
        self.mergerHalos = fits.BinTableHDU.from_columns(totalColumns).data
        
        return 1



    def calculateBulletVectors( self ):

        '''                                                                         
        Add to the fits structure, 'sources', the SI, SG and DI                     
        vectors, and also the beta values                                           
                                                                                
        This is taken from trajectory.pro used for my Science paper.                
        
        '''
        data = self.mergerHalos
        
        self.vector_sg = np.array([data['xGas'] - data['xStellar'], \
                          data['yGas'] -data['yStellar']])


        self.vector_sd = np.array([data['xDarkMatter'] -data['xStellar'], \
                          data['yDarkMatter'] -data['yStellar']])

        self.vector_gd =  np.array([data['xDarkMatter'] -data['xGas'], \
                          data['yDarkMatter'] -data['yGas']])

        self.vector_gs =  np.array([data['xStellar'] -data['xGas'], \
                          data['yStellar'] -data['yGas']])

        self.dist_sg = np.sqrt(np.sum( self.vector_sg**2, axis=0))
        self.dist_gi= np.sum(self.vector_gs*self.vector_gd, axis=0)\
          /np.sqrt(np.sum(self.vector_gs**2, axis=0))           # GS.GL/|GS|                                                  
        self.dist_si=np.sum(self.vector_sg*self.vector_sd, axis=0)/\
          np.sqrt(np.sum(self.vector_sg**2, axis=0))           # SG.SL/|SG|                                                   
        self.dist_di = (self.vector_sg[0]*self.vector_sd[1]-\
                            self.vector_sg[1]*self.vector_sd[0])\
          /np.sqrt(np.sum(self.vector_sg**2, axis=0))


        
        self.beta = self.dist_si / self.dist_sg

        self.betaPerp = self.dist_di / self.dist_sg
        


            
    def RemoveUnmatchedHalos( self ):

        if len(self.GasDarkMatter) > len(self.GasStellar):
            PrimaryArray = self.GasDarkMatter
            SecondaryArray = self.GasStellar
        else:
            PrimaryArray = self.GasStellar
            SecondaryArray = self.GasDarkMatter


        NewPrimaryArray =  \
          np.zeros(len(SecondaryArray) , dtype=PrimaryArray.dtype)
          
        SecondaryArrayIDs = \
          np.array([ i for i in SecondaryArray['ID_1']])

        NewCount = 0
        for iCount, iMergerPair in enumerate(PrimaryArray):
          
                        
            if (iMergerPair['ID_1'] in SecondaryArrayIDs):
                for iName in NewPrimaryArray.dtype.names:
                    NewPrimaryArray[NewCount][iName] = \
                      PrimaryArray[iCount][iName]
                NewCount += 1
                
        if len(self.GasDarkMatter) > len(self.GasStellar):
            self.GasDarkMatter = NewPrimaryArray
        else:
            self.GasStellar = NewPrimaryArray

                    
                    

    def CheckDarkMatterStellarDistance( self, distCut=10. ):

        
        xDarkMatterStellar = \
          self.GasDarkMatter['X_2'] - self.GasStellar['X_2']
        yDarkMatterStellar = \
          self.GasDarkMatter['Y_2'] - self.GasStellar['Y_2']
        dtype = self.GasDarkMatter.dtype
            
        distDarkMatterStellar = \
          np.sqrt( (xDarkMatterStellar)**2 + \
                       (yDarkMatterStellar)**2)

        nValid = len(distDarkMatterStellar[distDarkMatterStellar<distCut])
                                               
        NewGasStellar =np.zeros(nValid , dtype=dtype)

        NewGasDarkMatter = np.zeros(nValid , dtype=dtype)
        
        iValid = 0
        for iCount, iMergerPair in enumerate(self.GasDarkMatter):
            if distDarkMatterStellar[iCount] < distCut:
                
                for iName in dtype.names:
                    NewGasStellar[iValid][iName] = \
                      self.GasStellar[iCount][iName]
                      
                    NewGasDarkMatter[iValid][iName] = \
                      self.GasDarkMatter[iCount][iName]
                      
                iValid += 1
        
        if iValid ==0:
            return -4
        else:
            self.GasDarkMatter = NewGasDarkMatter
            self.GasStellar = NewGasStellar
            return 1
