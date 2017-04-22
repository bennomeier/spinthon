import itertools
import numpy as np
#import copy

from zeeman import zeemanProductBasis

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

        for k in range(dim/4, 2*dim/4):
            bFix = self.basis[k]
            self.basis[k] = 1/np.sqrt(2)*(self.basis[k + dim/4] + self.basis[k])
            self.basis[k + dim/4] = 1/np.sqrt(2)*(self.basis[k+dim/4] - bFix)

        #now lets rearange:

        self.basis = self.basis[2*dim/4:3*dim/4] + self.basis[0:2*dim/4] + self.basis[3*dim/4:]
            
        for b in self.basis:
            print b
