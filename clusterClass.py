import component_extractor as ce
import pyfits as fits
import numpy as np
import os as os
import matchMassComponents as mmc
import ipdb as pdb
from matplotlib import pyplot as plt

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
          '/Users/DavidHarvey/Documents/Work/Mergers/sims/BAHAMAS/KetelMount/BAHAMAS'
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

        if not os.path.isfile(self.GalaxyCat):
            pdb.set_trace()
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

        #self.UpdateStellarPositions()
        
        flag = self.checkCluster()
        if  flag < 0:
            return flag
        
        

        self.DarkMatterStars = \
          mmc.matchMassComponents( self.xDarkMatterPositions, \
                                   self.yDarkMatterPositions, \
                                   self.xStellarPositions, \
                                   self.yStellarPositions, \
                                   searchRadKpc = 40.)

        self.DarkMatterGas = \
          mmc.matchMassComponents( self.xDarkMatterPositions, \
                                   self.yDarkMatterPositions, \
                                    self.xGasPositions, \
                                   self.yGasPositions, \
                                   searchRadKpc = 200.)
        
        self.getGasStellar()

        if len(self.DarkMatterGas) == 0 & len(self.GasStellar) == 0:
            return -3                       
        if len(self.DarkMatterGas) != len(self.GasStellar):
            self.RemoveUnmatchedHalos()
        if len(self.DarkMatterGas) != len(self.GasStellar):
            raise ValueError("the remove unmatched didnt work and ive fucked up")
        distSDFlag = self.CheckDarkMatterStellarDistance()
        if distSDFlag < 0:
            return distSDFlag
        
        
        
        
        

        xGasColumn = \
          fits.Column( name='xGas', \
                           array=self.DarkMatterGas['X_2'], \
                           format='D')
                            
        yGasColumn = \
          fits.Column(  name='yGas', \
                            array=self.DarkMatterGas['Y_2'], \
                            format='D')
                            
        xDarkMatterColumn = \
               fits.Column(name='xDarkMatter', \
                            array=self.DarkMatterGas['X_1'], \
                            format='D')
        yDarkMatterColumn = \
               fits.Column(name='yDarkMatter', \
                           array=self.DarkMatterGas['Y_1'], \
                           format='D')

        xStellarColumn = \
          fits.Column( name='xStellar', \
                        array=self.GasStellar['X_2'], \
                        format='D')
                        
        yStellarColumn = \
          fits.Column(name='yStellar', \
                       array=self.GasStellar['Y_2'], \
                       format='D') 
        
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

        if len(self.DarkMatterGas) > len(self.GasStellar):
            PrimaryArray = self.DarkMatterGas
            SecondaryArray = self.GasStellar
        else:
            PrimaryArray = self.GasStellar
            SecondaryArray = self.DarkMatterGas


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
                
        if len(self.DarkMatterGas) > len(self.GasStellar):
            self.DarkMatterGas = NewPrimaryArray
        else:
            self.GasStellar = NewPrimaryArray

                    
                    

    def CheckDarkMatterStellarDistance( self, distCut=100. ):

        
        xDarkMatterStellar = \
          self.DarkMatterGas['X_2'] - self.GasStellar['X_2']
        yDarkMatterStellar = \
          self.DarkMatterGas['Y_2'] - self.GasStellar['Y_2']
        dtype = self.DarkMatterGas.dtype
            
        distDarkMatterStellar = \
          np.sqrt( (xDarkMatterStellar)**2 + \
                       (yDarkMatterStellar)**2)

        nValid = len(distDarkMatterStellar[distDarkMatterStellar<distCut])
                                               
        NewGasStellar =np.zeros(nValid , dtype=dtype)

        NewDarkMatterGas = np.zeros(nValid , dtype=dtype)
        
        iValid = 0
        for iCount, iMergerPair in enumerate(self.DarkMatterGas):
            if distDarkMatterStellar[iCount] < distCut:
                
                for iName in dtype.names:
                    NewGasStellar[iValid][iName] = \
                      self.GasStellar[iCount][iName]
                      
                    NewDarkMatterGas[iValid][iName] = \
                      self.DarkMatterGas[iCount][iName]
                      
                iValid += 1
        
        if iValid ==0:
            return -4
        else:
            self.DarkMatterGas = NewDarkMatterGas
            self.GasStellar = NewGasStellar
            return 1


    def getGasStellar( self ):
        '''
        Take the dark matter and find the stellar that is closest to it
        '''
        IdDarkMatterStars = self.DarkMatterStars['ID_1']

        xStellar = np.array([])
        yStellar = np.array([])
        IDStellar = np.array([])
        xGas = np.array([])
        yGas = np.array([])
        IDGas = np.array([])
        separation = np.array([])

        for iHalo in xrange(len(self.DarkMatterGas)):
            index = IdDarkMatterStars == self.DarkMatterGas['ID_1'][iHalo]
            
            xStellar = np.append(xStellar,  \
              self.DarkMatterStars['X_2'][ index ])
            yStellar =  np.append(yStellar, \
              self.DarkMatterStars['Y_2'][ index ])
            IDStellar =  np.append(IDStellar, \
              self.DarkMatterStars['ID_2'][ index ])
              
            iSeparation = np.sqrt( (self.DarkMatterStars['X_2'][ index ] - \
                            self.DarkMatterGas['X_2'][iHalo])**2 + \
                       (self.DarkMatterStars['Y_2'][ index ] - \
                            self.DarkMatterGas['Y_2'][iHalo])**2   )
            
                            
            separation = np.append( separation, iSeparation)
              

            xGas = np.append( xGas, self.DarkMatterGas['X_2'][iHalo])
            yGas = np.append( yGas, self.DarkMatterGas['Y_2'][iHalo])
            IDGas = np.append( IDGas, self.DarkMatterGas['ID_2'][iHalo])
            

        X_1 = fits.Column(name='X_1', array=xGas, format='D')
         
        Y_1 = fits.Column(name='Y_1', array=yGas, format='D')

        ID_1 = fits.Column(name='ID_1', array=IDGas, format='D')


        X_2 = fits.Column(name='X_2', array=xStellar, format='D')
         
        Y_2 = fits.Column(name='Y_2', array=yStellar, format='D')

        ID_2 = fits.Column(name='ID_2', array=IDStellar, \
                              format='D')

        separationCol =  fits.Column(name='Separation', array=separation, \
                              format='D')


        columns = [X_1, Y_1, ID_1, X_2, Y_2, ID_2, separationCol]
        
        self.GasStellar = fits.BinTableHDU.from_columns(columns).data
        if len(self.GasStellar) != len(self.DarkMatterGas):
            raise ValueError("Gas vector not the correct length")


    def UpdateStellarPositions( self, massCut=0. ):
        self.getClusterMembers()
        massCutIndex = self.clusterMembers['mass100'] > massCut
        self.xStellarPositions = \
          self.clusterMembers['x'][ massCutIndex ]

        self.yStellarPositions = \
          self.clusterMembers['y'][ massCutIndex ]
    


    def plotPositions(self):

        plt.plot(self.xGasPositions, self.yGasPositions, 'r.')
        plt.plot(self.xDarkMatterPositions, self.yDarkMatterPositions, 'b.')
        plt.plot(self.xStellarPositions, self.yStellarPositions, 'g.')



    def BinaryCluster( self, pixelScale=5., aperature=200, closeCut=15,\
                        massContrastCut=0.5, massCut=11.5, GasDistCut=50.):
        '''
        Retur a boolean if there is a galaxy of mass contrast to the BCG
        greater than massContrastCut (default 0.5) and absolute mass 
        of greater than 11.5
        within aperature (default =  100kpc ) of the central BCG 
        (projected distance)
        and not wihtin the closeCUt (the prior of the BCG)
            
        '''
        self.getClusterMembers()
        massCutIndex = self.clusterMembers['mass100'] > massCut
        BinaryHalo = []
        galaxyMass = self.clusterMembers['mass100']
        for iHalo in self.mergerHalos:

            distance = np.sqrt( (iHalo['xStellar'] - self.clusterMembers['x'])**2 + \
                            (iHalo['yStellar'] - self.clusterMembers['y'])**2)*pixelScale

            nearestBCGmass = self.clusterMembers['mass100'][ np.argmin(distance) ]
        
            MassContrast = 10**nearestBCGmass / \
                                10**self.clusterMembers['mass100']



            DistanceBool = (distance < aperature) &  (distance > 0)
            MassContrastBool = (MassContrast[DistanceBool] > massContrastCut)
            MassBool =  (galaxyMass[DistanceBool] > massCut)
            CloseBool = (distance[ DistanceBool ] < closeCut)
    
            BinaryHalo.append(np.any(MassContrastBool | MassBool | CloseBool))

        '''
        BinaryHalo = []
        for iHalo in self.mergerHalos:

            distance = np.sqrt( (iHalo['xDarkMatter'] - self.mergerHalos['xDarkMatter'])**2 + \
                            (iHalo['yDarkMatter'] - self.mergerHalos['yDarkMatter'])**2)

            BinaryHalo.append(np.any( (distance > 0.) & (distance < GasDistCut)))
        '''
        return np.array(BinaryHalo) == False
