import itertools

import numpy as np
import spinthon.spin.single.data
import spinthon.spin.single.operators

from ..basis import zeeman
from ..basis.vector import ElementaryKet, ElementaryBra


class Spin(object):
    """represent a single spin"""
    def __init__(self, name):
        self.gamma = spinthon.spin.single.data.gammaListAll[name]
        self.spin = spinthon.spin.single.data.spinListAll[name]

        """this returns the matrix representation of a single spin 
        in a single spin hilbert state with elementary basis kets."""
        self.spinOps = spinthon.spin.single.operators.getSpinOperators(self.spin)

    def __getattr__(self, attr):
        #allow retrieval of e.g. Ix by calling spin.Ix
        return self.spinOps[attr]


class spinSystem(object):
    def __init__(self, spins, verbose = False, basis = "zeeman"):
        """Spins: A list of spins, e.g. ["1H", "1H", "1H", "13C"]"""

        self.spinSystem = []

        self.dimension = 1
        
        for s in spins:
            S = Spin(s)
            self.spinSystem.append(S)

            self.dimension = int(self.dimension*(2*S.spin+1))

        maxVal = tuple([s.spin for s in self.spinSystem])
        ranges = [np.arange(-s.spin, s.spin + 1) for s in self.spinSystem]

        elementaryStatesList = list(itertools.product(*ranges))
            
        self.ekets = [ElementaryKet(l, maxVal) for l in elementaryStatesList]
        self.ebras = [ElementaryBra(l, maxVal) for l in elementaryStatesList]
            
        if basis == "zeeman":
            self.basis = zeeman.zeemanProductBasis(self)
            
        if verbose:
            print("Spin System initialized.")
            print("Dimensionality: {}".format(self.dimension))

            
        
    def __getitem__(self, pos):
        """support indexing of the spin system"""
        return self.spinSystem[pos]


        
if __name__ == "__main__":

    methyl = spinSystem(["1H", "1H", "1H"])
    
    methyl13C = spinSystem(["1H", "1H", "1H", "13C"])

    h2o17 = spinSystem(["1H", "1H", "17O"])
    
