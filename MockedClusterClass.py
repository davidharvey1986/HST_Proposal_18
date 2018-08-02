import numpy as np
from matplotlib import pyplot as plt
import matchMassComponents as mmc
import pyfits as fits

class MockedCluster:

    def __init__( self, \
                    DarkMatterStellarOffset,  \
                    GasStellarOffset,
                    DarkMatterError, \
                    GasError, \
                    nSpuriousHalos):

        '''
        Create a mocked cluster sample
        DarkMatterStellarOffsert = offset between dm and stars
        '''
        
        self.nSpuriousHalos  = nSpuriousHalos
        self.DarkMatterStellarOffset = \
          DarkMatterStellarOffset
        self.GasStellarOffset = \
              GasStellarOffset
        self.GasError = \
              GasError
        self.DarkMatterError = \
              DarkMatterError
        self.xClusterCenter = 1000.
        self.yClusterCenter = 1000.
        self.nHalos = 2
        self.getHaloPosition()
        self.getGasPositions()
        self.getDarkMatterPositions()
        self.getStellarPositions()
        self.rotateSystem( )
        self.addNoise( )
        self.getClusterMass()

    def addNoise( self ):
        '''
        Add Random halos to the positions

        '''

        #self.xGasPositions = np.append(self.xGasPositions, \
          #                                 self.RandomPositions())
        #self.yGasPositions = np.append(self.yGasPositions, \
         #                                  self.RandomPositions())
        self.xDarkMatterPositions = \
          np.append(self.xDarkMatterPositions,  self.RandomPositions())
        self.yDarkMatterPositions = \
          np.append(self.yDarkMatterPositions, self.RandomPositions( ))
        self.xStellarPositions = np.append(self.xStellarPositions, \
                                               self.RandomPositions( ))
        self.yStellarPositions = np.append(self.yStellarPositions, \
                                               self.RandomPositions( ))

    def RandomPositions(self ):

        return  np.random.uniform(0.,2000., self.nSpuriousHalos)

        
    def rotateSystem( self ):

        
        self.RotationAngle = np.random.uniform( 0, 2.*np.pi, self.nHalos )
                     
        self.xGasPositions, self.yGasPositions = \
          self.rotationMatrix( self.xGasPositions, self.yGasPositions)

        self.xDarkMatterPositions, self.yDarkMatterPositions = \
          self.rotationMatrix( self.xDarkMatterPositions, \
                               self.yDarkMatterPositions)

        
    def rotationMatrix( self, x, y ):

        diffx = x - self.xStellarPositions
        diffy = y - self.yStellarPositions
        
        NewX = diffx*np.cos(self.RotationAngle) - \
          diffy*np.sin(self.RotationAngle)

        NewY = diffx*np.sin(self.RotationAngle) + \
          diffy*np.cos(self.RotationAngle)
          
        
        NewX += self.xStellarPositions
        NewY += self.yStellarPositions

        return NewX, NewY
    
    def getHaloPosition( self ):
        self.xHaloPosition = np.random.uniform(0, 2000, self.nHalos)
        self.yHaloPosition = np.random.uniform(0, 2000, self.nHalos)

    def getGasPositions( self ):

        self.xGasPositions = self.GasStellarOffset + \
            np.random.normal(self.xHaloPosition, \
                             self.GasError, \
                             self.nHalos) 

        self.yGasPositions = \
          np.random.normal(self.yHaloPosition, \
                           self.GasError, \
                           self.nHalos) 
          

    def getDarkMatterPositions( self ):
        self.xDarkMatterPositions = self.DarkMatterStellarOffset + \
            np.random.normal(self.xHaloPosition , \
                             self.DarkMatterError, \
                             self.nHalos) 
            
        self.yDarkMatterPositions = \
            np.random.normal(self.yHaloPosition , \
                             self.DarkMatterError,  \
                             self.nHalos) 
            
          
    def getStellarPositions( self ):

        self.xStellarPositions = self.xHaloPosition
        self.yStellarPositions = self.yHaloPosition

    def getClusterMass( self ):

        self.mass = 14.6

          
    def combineMassComponents( self ):
        #Combine the three mass components
        
        self.GasDarkMatter = \
            mmc.matchMassComponents( self.xGasPositions, \
                                   self.yGasPositions, \
                                   self.xDarkMatterPositions, \
                                   self.yDarkMatterPositions, \
                                   searchRadKpc = 100.)

        self.GasStellar = \
            mmc.matchMassComponents( self.xGasPositions, \
                                   self.yGasPositions, \
                                   self.xStellarPositions, \
                                   self.yStellarPositions, \
                                   searchRadKpc = 100.)



        
        
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
                           format='D')
                            
        yGasColumn = \
          fits.Column(  name='yGas', \
                            array=self.GasDarkMatter['Y_1'], \
                            format='D')
                            
        xDarkMatterColumn = \
          fits.Column(name='xDarkMatter', \
                            array=self.GasDarkMatter['X_2'], \
                            format='D')
        yDarkMatterColumn = \
          fits.Column(name='yDarkMatter', \
                           array=self.GasDarkMatter['Y_2'], \
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

                    
                    

    def CheckDarkMatterStellarDistance( self, distCut=20. ):
       
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


      

          

          
      
