import itertools
import numpy as np

import vector
reload(vector)

from vector import BasisKet, Ket



class zeemanProductBasis(object):
    def __init__(self, spinSystem):
        maxVal = tuple([s.spin for s in spinSystem])
        ranges = [np.arange(-s.spin, s.spin + 1) for s in spinSystem]
        
        self.name = "Zeeman Product Basis"
        self.basisketsList = list(itertools.product(*ranges))
        
        self.basiskets = [BasisKet(l, maxVal) for l in self.basisketsList]

        #for the zeeman basis the coefficients are just one in one position and zeros elsewhere

        self.basis = []

        for k in range(spinSystem.dimension):
            coefficients = np.zeros(spinSystem.dimension)
            coefficients[k] = 1

            #note that it maybe slightly redundant to store the basiskets.
            self.basis.append(Ket(coefficients, self.basiskets))

        
        
    def __getitem__(self, pos):
        """support indexing of the basis"""
        return self.basis[pos]
        
