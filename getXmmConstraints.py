'''
The idea is to use the mergers as a way to constrain the sterile
neutrino through the 3.5 keV line.

To do this i need to estimate the constraints

What is the esimated flux of these clusters?
What constraints on the 3.5keV line will we get?
How many offsets between dark matter and gas do we expect?

'''
import ciao_contrib.runtool as ciao

import pyfits as fits
from HSTconstraints import *
import idlsave as idl
from astro_tools_extra import ra_separation as raSep
from lensing.lensing_parameters import *
import pyRRG as pyRRG
import RRGtools as at
import pysex as sex

def getDistributionStarsGas():
    '''
    Distrubtion in stars and gas
    '''

    harveyDataDir = \
      '/Users/DavidHarvey/Documents/Work/SIDM_science/clusters/'

    clusters = glob.glob(harveyDataDir+'/*')
    listOfMergerHalos = []
    for iCluster in clusters:
        print iCluster
        if os.path.isfile( iCluster+'/bulletsFLUX.sav'):
            
            halos = idl.read(iCluster+'/bulletsFLUX.sav')['halos']
            for iHalo in halos:
                clusterName = iCluster.split('/')[-1]
                iMerger = mergerHalo(clusterName, iHalo)
                iMerger.calculateBulletVectors()
                listOfMergerHalos.append(iMerger)



    simName = 'CDM'
    #simulatedClusters = ClusterSample(simName)

    #simulatedClusters.extractMergerHalos()

    #simulatedClusters.CalculateOffsetVectors(nClusters=5)

    bins = np.linspace(0,1000,100)
    
    totalSGdistribution = [ i.dist_sg[0] for i in listOfMergerHalos]
    
    proposedClusterRadius = getProposedClusterRadii()
    y,x,p = plt.hist( totalSGdistribution, bins=bins, density=True)
    plt.hist( proposedClusterRadius, bins=bins,density=True)
    plt.plot( [np.median(proposedClusterRadius),np.median(proposedClusterRadius)],[0.,np.max(y)],'-')

    plt.show()

def getProposedClusterRadii():
    '''
    ge the radius of all the clsuters that i have
    proposed (from wavdetect file, so not ideal
    '''
    pklFile = "Pickles/proposedClusterRadii.pkl"
    if os.path.isfile(pklFile):
        radii = np.array(pkl.load(open(pklFile,'rb')))
        return radii[radii != 100]
    
    
    #Get the mean radius of each halo in the same
    dtype = [('Name',object), ('Exp', int),('Redshift',float),\
            ('Roll',float), ('vmax', float), ('nHalos',int),('Xmm',float)]

    
    proposedClustersDir = \
      '/Users/DavidHarvey/Documents/Papers/Hubble_Proposal/proposal_17_resub_harvey/data'
    proposedClusters = \
      np.loadtxt( proposedClustersDir+'/allClusters.lis', dtype=dtype)
    haloSize=[]
    for iCluster, iClusterName in enumerate(proposedClusters['Name']):
        sourceFile = \
          glob.glob(proposedClustersDir+'/'+iClusterName+'*/Wavdetect/extras/'+\
                      iClusterName+'*_source.fits')[0]

        clusterSources = fits.open(sourceFile)[1].data
        wavFile = \
          glob.glob(proposedClustersDir+'/'+iClusterName+'*/Wavdetect/'+\
                      iClusterName+'*_wavdetect.fits')[0]

        pixelSizeArcsec = fits.open(wavFile)[0].header['CDELT2'] *3600.
        pixelSizeKpc = pixelSizeArcsec/206265.* \
          (ang_distance( proposedClusters['Redshift'][iCluster] )*1e3)
        #sort htem in significance
        argSortedHalos = \
          np.argsort( clusterSources['SRC_SIGNIFICANCE'] )[::-1]
        clusterSources=  clusterSources[argSortedHalos]
        nHalos = \
          np.min([ proposedClusters['nHalos'][iCluster], \
                len(clusterSources['R'])])
        for iSource in xrange(nHalos):
            sizeKpc =  clusterSizeSex( wavFile, clusterSources[iSource])
            haloSize.append(sizeKpc)
    pkl.dump(haloSize, open(pklFile, 'wb'))
    return haloSize



    
class mergerHalo:
    def __init__( self, clusterName, halo ):
        self.clusterName = clusterName
        self.mergerHalo = { \
            'xGas':halo['GAS']['RA'], \
            'yGas':halo['GAS']['DEC'], \
            'xStellar':halo['GAL']['RA'], \
            'yStellar':halo['GAL']['DEC'], \
            'xDarkMatter':halo['DM']['RA'], \
            'yDarkMatter':halo['DM']['DEC'] }
            
        self.getRedshift()
        self.getAngularDistance()

    def getAngularDistance( self ):
        #ang distance in kpc
        self.angularDistance = \
          ang_distance(self.redshift)*1e3

        #factor to convert arcseconds to kpc
        self.arcsec2kpc = \
          1./206265*self.angularDistance
          
        
    def getRedshift( self ):
        
        dataDir = '/Users/DavidHarvey/Documents/Work/SIDM_science/'
        dtype = [('name',object), ('info', object), ('z_lens', float)]
        redshiftInfo = \
          np.loadtxt(dataDir+'/data/redshiftInformation.lis', \
                         dtype=dtype)
        clusterInfoNames = \
          np.array([ i.split('_')[0] for i in redshiftInfo['name']])
        index = self.clusterName == clusterInfoNames
        self.redshift = redshiftInfo['z_lens'][ index ][0]
        
        
    def calculateBulletVectors( self ):

        '''                                                                         
        Add to the fits structure, 'sources', the SI, SG and DI                     
        vectors, and also the beta values                                           
                                                                                
        This is taken from trajectory.pro used for my Science paper.                
        
        '''
        data = self.mergerHalo
        
        self.vector_sg = \
          np.array([raSep(data['xGas'], data['yGas'],\
                          data['xStellar'], data['yGas']), \
                          (data['yGas'] -data['yStellar'])*3600.])\
                          *self.arcsec2kpc


        self.vector_sd = \
          np.array([raSep(data['xDarkMatter'], data['yDarkMatter'],\
                        data['xStellar'],data['yDarkMatter']), \
                          (data['yDarkMatter'] -data['yStellar'])*3600.])\
                          *self.arcsec2kpc

        self.vector_gd =  \
          np.array([raSep(data['xDarkMatter'], data['yDarkMatter'],\
                          data['xGas'], data['yDarkMatter']), \
                          (data['yDarkMatter'] -data['yGas'])*3600.])\
                          *self.arcsec2kpc

        self.vector_gs =  \
          np.array([raSep(data['xStellar'],data['yStellar'],\
                          data['xGas'], data['yStellar']), \
                          (data['yStellar'] -data['yGas'])*3600.])\
                          *self.arcsec2kpc

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
        
def xray_pix2sky( inFile, x, y, asol='none'):
    ra = []
    dec = []
    ciao.dmcoords.punlearn()
    for i in xrange(len(x)):
	dm = ciao.dmcoords( infile=inFile, asol=asol, \
                            x=x[i], y=y[i], opt='sky', celfmt='deg')


        ra.append( np.float(ciao.dmcoords.ra))
	dec.append( np.float(ciao.dmcoords.dec))

    return np.array(ra), np.array(dec)


def clusterSizeSex( image, clusterSource):
    '''
    Get cluster size in arcseconds
    '''
    try:
        sources = sex.run(image)

        returnSource = sources[ np.argmin(raSep(sources['X_WORLD'],sources['Y_WORLD'],clusterSource['RA'],clusterSource['DEC'],abs=True))]
                                          
        radius = returnSource['FLUX_RADIUS1']*fits.open(image)[0].header['CDELT2']*3600/206265*ang_distance(0.2)*1e3
    except:
        radius =100

    return radius
                                   
def clusterSizeWavdetect( image, clusterSource):
    '''
    Get cluster size in arcseconds
    '''
    sourceSize = np.sqrt(clusterSource['R'][0]*clusterSource['R'][1])
    xImagePoints = np.array([ clusterSource['X'], clusterSource['X']])
    yImagePoints = np.array([ clusterSource['Y'], clusterSource['Y']+sourceSize])

    pointRA, pointDEC =  xray_pix2sky( image, xImagePoints, yImagePoints )
    size = np.abs(clusterSource['DEC'] - pointDEC[1])*3600
    return size

