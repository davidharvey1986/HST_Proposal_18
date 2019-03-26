from scipy import optimize
from matplotlib import pyplot as plt
import numpy as np
import ipdb as pdb

def plotTrend( x, y, yError, simName, plotLims=None, color=None):

    #IF i want to fit all together, and regularise (very complex)                                            
    #then:

   
    trendParams, pError = getTrend( x, y, yError )

   
    print y
    ax = plt.gca()
    if plotLims is None:
        plotLims = [np.min(x), np.max(x)]
    xPlot = np.linspace(plotLims[0], plotLims[1], 100)
    ax.plot( xPlot, fitfunc( xPlot, *trendParams ), '--', color=color)
    
    pLower = [trendParams[0]+pError[0], \
                  trendParams[1]-pError[1]]
    pUpper = [trendParams[0]-pError[0], \
                  trendParams[1]+pError[1]]
    

    ax.fill_between( xPlot, fitfunc( xPlot, *pUpper), \
                         fitfunc( xPlot, *pLower), \
                         alpha=0.3, color=color )

    return trendParams
def fitfunc( x, p1, p2):
    return p1 + p2*x
def getTrend( x, y, yError ):

    '''                                                                                                      
    Fit a straight line to the data                                                                          
    and return with an error                                                                                 
    '''

    pinit = [1.0, -1.0]
    args = ( x, y, yError )
    pfinal, covar = \
      optimize.curve_fit(fitfunc, x, y, sigma=yError)

    print("y = %0.4f + %0.4f * x " % tuple(pfinal))
    pError = [ np.sqrt( covar[0][0] ), np.sqrt( covar[1][1] ) ]

    return pfinal, pError
