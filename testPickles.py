'''
A lot of the pickles seem fucked so need to test them
'''
import glob
import pickle as pkl
from clusterComponents import *

def main():
    pickleFiles = glob.glob('Pickles/StellarAndDMPositions/*')
    i=0
    for iPkl in pickleFiles:
        pklObject = pkl.load(open(iPkl,'rb'))
        try:
            pklObject = pkl.load(open(iPkl,'rb'))
        except:
            i+=1
            print iPkl

    
    print i,len(pickleFiles)
