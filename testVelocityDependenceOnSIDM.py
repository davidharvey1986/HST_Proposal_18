
import lensing_parameters as l
import emcee
import fitVelocityDependentCrossSection as fitVDcross
import numpy as np
import corner 
from matplotlib import pyplot as plt
def main( ):
    '''
    Can we constrain velocity
    '''
    #Set up the observable
    nClusters = 130
    nMassBins = 5
    nClustersPerBin = nClusters / nMassBins

    
    #Theses are the clusters i will observe
    mass = np.linspace(14,15.5,nMassBins)
    redshift = 0.2
    rad = l.virial_radius( 10**mass, redshift)*1e3
    velocity = np.sqrt(10**mass*l.grav_const()/rad)
    print(velocity)
    #shear goes in and sis goes as propto einstein rad =priot m^0.5
    #so signal to noise will go as propto velocity
    signalToNoise = velocity / np.max(velocity )

    
    errorOnCrossSectionPerBin = signalToNoise / np.sqrt(nClustersPerBin)
    
    signal = fitVDcross.getVelocityDependentCrossSection( velocity, 3.04, 560.)*0.



    inputDict = {'x':velocity, 'y':signal, 'error':errorOnCrossSectionPerBin}
    
    fitClass = fitVDcross.fitVelocityDependentCrossSection( inputDict )
    corner.corner( fitClass.samples, range=[[0.,5.],[0.,1000]])
    plt.show()

if __name__ == '__main__':
    main()
    


