import pyfits as fits
import numpy as np
import os as os

def matchComponents( self, xCatA, yCatA, \
                             xCatB, yCatB, \
                             searchRadKpc=1. ):
    '''
    Use the classic matching program
    but to do this i will create fake RA and DEC
    '''

        
    catAFitsTable = \
      self.createPosFitsTable( self.xCatA, \
                                    self.yCatA)
                                   

    catBFitsTable = \
      self.createPosFitsTable( self.xCatB, \
                                   self.yCatB)

    randomStr = str(np.random.random_integers(0,10000))
    gasFitsTable.writeto('catA_'+randomStr+'.fits')
    dmFitsTable.writeto('catB_'+randomStr+'.fits')
    
    matchedCat = \
      self.runMatch( 'catA_'+randomStr+'.fits', \
                                'catB_'+randomStr+'.fits', \
                                 search_rad=searchRadKpc)[1].data
       

    os.system('rm -fr catA_'+randomStr+'.fits')
    os.system('rm -fr catB_'+randomStr+'.fits')
    
    return matchedCat
    
def createPosFitsTable( self, x, y ):

    #need to normallise x and y to ra and dec limits
        
    RA = fits.Column(name='X', array=x, format='K')
                             
    DEC = fits.Column(name='Y', array=y, format='K')

    ID = fits.Column(name='ID', array=np.arange(len(x)), \
                              format='K')
                              
    
    return fits.BinTableHDU.from_columns([RA,DEC,ID])
    


def runMatch( self, catA, catB, search_rad=20. ):
    stilts_path='/Users/DavidHarvey/Library/Stilts' 
        
    command_str =stilts_path+'/stilts.sh tmatch2 '\
          +'in1="'+catA+'" in2="'+catB+'" '\
          +'matcher=2d values1="X Y" values2="X Y" '\
          +'params="'+str(search_rad)+'" '\
          +'out=matched_A_B.fits'


    os.system(command_str)

    matched_cat = fits.open('matched_A_B.fits')

    os.system("rm -fr matched_A_B.fits")

    return matched_cat