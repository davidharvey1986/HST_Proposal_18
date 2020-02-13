


def velocityDependentCrossSection( velocity, velocityBreak, sigma0=3.):
    
   return sigma0/(1+ (velocity/velocityBreak)**2)**2

