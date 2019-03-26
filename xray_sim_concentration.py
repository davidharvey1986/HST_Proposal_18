import pyfits as fits
import numpy as np
import ipdb as pdb
from matplotlib import pyplot as plt
import lensing as lens
from scipy.ndimage.filters import gaussian_filter as gauss
def xray_sim_concentration( cluster_number, r_grid, \
                                sim = 'AGN_TUNED_nu0_L400N1024_WMAP9',\
                                clusterRedshiftStr='z_0.250'):
    '''
    Get the ratio of the emiision within 100 kpc / 400 kpc

    r_grid : a array nxn dimensions where n is the size of the image
    file, (2000x2000) with each pixel the radial distance from the centre.

   '''
   
        

    
    dataDir = '/Users/DavidHarvey/Documents/Work/Mergers/sims/BAHAMAS/KetelMount/BAHAMAS/'+sim+'/'+clusterRedshiftStr
    cluster_catalogue = np.loadtxt( dataDir+'/catalog.dat',
                                    dtype=[('id', int), ('fof', int), ('mass', float), ('r200', float)])
    
    cluster_mass = 10**(cluster_catalogue['mass'][ cluster_catalogue['id'] == cluster_number])
    cluster_file = dataDir+'/MAPS/cluster_'+str(cluster_number)+'_xray.fits'

    xray_data = fits.open( cluster_file )[0].data



    return np.sum(xray_data[ r_grid < 100 ]) / np.sum(xray_data[ r_grid < 400])

def get_all_xray_concentrations( clusterInts, x0=1000, y0=1000, \
                                     sim='AGN_TUNED_nu0_L400N1024_WMAP9',\
                                     clusterRedshiftStr='z_0.250' ):
    '''
    clusterInts : the integers of each cluster name
    loop through all xray maps and get concentrations
    '''

    dataDir = '/Users/DavidHarvey/Documents/Work/Mergers/sims/BAHAMAS/KetelMount/BAHAMAS/'+sim+'/'+clusterRedshiftStr
    image_size = fits.open( dataDir+'/MAPS/cluster_'+clusterInts[0]+'_xray.fits' )[0].data.shape

    pixel_size = 5 #kpc
    x_vector = np.arange( image_size[0] )*pixel_size - x0*pixel_size
    y_vector = np.arange( image_size[1] )*pixel_size - y0*pixel_size
    x_grid, y_grid = np.meshgrid( x_vector, y_vector )
    r_grid = np.sqrt( x_grid**2 + y_grid**2)

    xray_concentration = []
    for i, iCl in enumerate(clusterInts):
        print("Cluster %s (%i/%i)" % (iCl,i+1,len(clusterInts)))
        iConcentration = xray_sim_concentration( iCl, r_grid, sim=sim,\
                                clusterRedshiftStr=clusterRedshiftStr)
        xray_concentration.append( iConcentration )

    return xray_concentration
    
