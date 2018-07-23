import pysex as sex
from scipy.ndimage.filters import gaussian_filter as gauss_filter
import pyfits as fits
from matplotlib import pyplot as plt
import numpy as np
from numpy.lib.recfunctions import append_fields as append_rec
import lensing as lensing
import ipdb as pdb
import os as os

def component_extractor( image_file, smoothing=20, \
                             pixel_size_kpc=5, redshift=0.25,\
                             filtername='gauss_5.0_9x9.conv'):
    '''
    This is a general script that should extract a given
    mass component from the simulations.   

    
    '''

    image = fits.open( image_file )[0].data

    if smoothing is not None:
        smoothed_image = gauss_filter( image, smoothing )
    else:
        smoothed_image = image
    pixel_size_degrees = pixel_size_kpc/(lensing.ang_distance( redshift )*1e3)*180./np.pi

    #Create a random string to component find
    file_name = 'component_finding_'+str(np.random.random_integers(0,1000))+'.fits'
    fits.writeto( file_name, smoothed_image, clobber=True )
    
    conf_path = '/Users/DavidHarvey/Library/Code/IDL/rrg/code/shape/sex_files/'
    hot_conf = conf_path+'wavdetectpy.sex'

    filtername = conf_path+'/'+filtername
   
    conf_args = {'WEIGHT_TYPE':'NONE',
                 'CLEAN':'Y',
                 'DEBLEND_NTHRESH':4, 
                 'FILTER_NAME':filtername,
                    'FILTER':'YES'}
    
    sources = sex.run(  file_name,   conf_file=hot_conf, conf_args=conf_args )
   
    os.system('rm -fr '+file_name)
    ra = sources['X_IMAGE']*pixel_size_degrees
    dec = sources['Y_IMAGE']*pixel_size_degrees
    sources = append_rec( sources, 'RA', ra, usemask=False, asrecarray=True)
    sources = append_rec( sources, 'DEC', dec)

    
    
    return sources

    
