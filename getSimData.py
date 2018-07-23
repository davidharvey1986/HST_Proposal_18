'''

A single script to get the sim data

its uses BCG_check_centre, but it has become all very
convoluted so this is a wrapper

'''


def getSimNameList():
    Sims = ['CDM','SIDM0.1','SIDM', 'CDMhires','SIDMhires','obs']
    return Sims


      
def getSimList( SimName ):
    if SimName == 'CDM':
        Sims = [ 'AGN_TUNED_nu0_L400N1024_WMAP9_LOMASS',
                    'AGN_TUNED_nu0_L400N1024_WMAP9' ]
    if SimName == 'CDMhires':
        Sims = [  'AGN_RECAL_nu0_L100N512_WMAP9']
    if SimName == 'SIDMhires':
        Sims = [ 'AGN_RECAL_nu0_SIDM_sig1_L100N512_WMAP9' ]
    if SimName == 'CDM':
        Sims = [       'AGN_TUNED_nu0_L400N1024_WMAP9_LOMASS',
                    'AGN_TUNED_nu0_L400N1024_WMAP9' ]
    if SimName == 'SIDM':
        Sims = [       'AGN_TUNED_nu0_SIDM_L400N1024_WMAP9_LOMASS',
                   'AGN_TUNED_nu0_SIDM_L400N1024_WMAP9' ]
    if SimName == 'SIDM0.1':
        Sims = [ 'AGN_TUNED_nu0_SIDM_sig0p1_L400N1024_WMAP9_LOMASS',\
                     'AGN_TUNED_nu0_SIDM_sig0p1_L400N1024_WMAP9']

    if SimName == 'obs':
        Sims = [ 'obs' ]

    return Sims


