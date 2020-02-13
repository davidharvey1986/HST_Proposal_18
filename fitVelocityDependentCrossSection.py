'''
This script fits a velocity dependent crosssection

'''
import pickle as pkl
from scipy.optimize import curve_fit
import emcee
from scipy.stats import norm
from scipy.stats import chi2

import numpy as np
import os
import ipdb as pdb
from matplotlib import pyplot as plt


def lnprob( theta, xTrue, yTrue, error ):
   
    yTheory = getVelocityDependentCrossSection( xTrue, *theta)


    prob = np.sum(norm.logpdf( yTheory, yTrue, scale=error))
    
    if np.isnan(prob):
        pdb.set_trace()
        return -np.inf
    if (theta[0] < 0) | (theta[0] > 10):
        return -np.inf
    if theta[1] < 0:
        return -np.inf
    return prob
    
   
def getVelocityDependentCrossSection( velocity, sigma0, velocityBreak):
    
        return sigma0/(1+ (velocity/velocityBreak)**2)**2

class fitVelocityDependentCrossSection:
    
    def __init__( self, pdf ):
        '''
        Init the pdf 
          
        the pdf is a dict of 'x', 'y', 'yLensPlane', 'yError','yLensPlaneError'
        
        '''                   
        self.pdf = pdf
       
        self.fitCrossSection()



    def fitCrossSection( self, nthreads=4  ):

  
        nwalkers = 20

        ndim = 2
        burn_len=100
        chain_len=1000
        pos0 = np.random.rand(nwalkers,ndim)
        pos0[:,1] = np.random.rand( nwalkers )*1000
        
        
        args = (self.pdf['x'], self.pdf['y'], \
                    self.pdf['error'] )
    
        dmsampler = emcee.EnsembleSampler(nwalkers, ndim, \
                                            lnprob, \
                                          args=args, \
                                          threads=nthreads)
                                          
        pos, prob, state  = dmsampler.run_mcmc(pos0, burn_len)


    
        pos, prob, state  = dmsampler.run_mcmc(pos, chain_len,\
                                        progress=True)
        self.samples = dmsampler.flatchain

        errorLower, median, errorUpper = \
          np.percentile(self.samples, [16, 50, 84], axis=0)

        error = np.mean([median - errorLower, errorUpper - median], axis=0)

        self.params = {'params':median, 'error':error}


   
    def getVelocityDependentCrossSection( velocity, velocityBreak, sigma0=3.):
    
        return sigma0/(1+ (velocity/velocityBreak)**2)**2



    def saveSamples( self, pklFile):
        pkl.dump( self.samples, open(pklFile, 'wb'))
