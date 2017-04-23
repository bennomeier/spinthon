import itertools
import numpy as np
#import copy

from .zeeman import zeemanProductBasis

class singletTripletBasis(zeemanProductBasis):
    def __init__(self, spinSystem):
        zeemanProductBasis.__init__(self, spinSystem)

        #assume the singlet-triplet quantum numbers come first.
        #hence part all the kets into four sets of the same length

        #the first ones are -1/2, -1/2,x 
        #the second ones are -1/2, 1/2, x
        #the third ones are 1/2, -1/2, x
        #the fourth ones are 1/2, 1/2, x

        self.name = "Singlet Triplet Basis"

        #hence first and fourths onces are triplets, second and third ones need to be mixed 

        dim = spinSystem.dimension

        start = int(dim/4)
        stop = int(2*dim/4)
        
        for k in range(start, stop):
            bFix = self.Kets[k]
            self.Kets[k] = 1/np.sqrt(2)*(self.Kets[k + start] + self.Kets[k])
            self.Kets[k + start] = 1/np.sqrt(2)*(self.Kets[k+start] - bFix)

        #now lets rearange:

        self.Kets = self.Kets[2*start:3*start] + self.Kets[0:2*start] + self.Kets[3*start:]
        self.Bras = [k.getBra() for k in self.Kets]
        
        for b in self.Bras:
            print(b)
