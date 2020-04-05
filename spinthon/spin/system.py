import itertools

import numpy as np
from scipy.linalg import expm

import spindata
import spinthon.spin.single.operators
import spinthon.spin.geometry

from ..basis import zeeman
from ..basis.vector import ElementaryKet, ElementaryBra


class Spin(object):
    """represent a single spin"""
    def __init__(self, name):
        self.gamma = spindata.gamma(name)
        self.spin = spindata.spin(name)

        """this returns the matrix representation of a single spin 
        in a single spin hilbert state with elementary basis kets."""
        self.spinOps = spinthon.spin.single.operators.getSpinOperators(self.spin)

    def rotationOp(self, beta, phi):
        """This returns the single spin  matrix representation of the rotation operator
        - beta: flip angle
        - phi: phase
        """
        return expm(-1j*beta*(self.spinOps["Ix"]*np.cos(phi) + self.spinOps["Iy"]*np.sin(phi)))
 
    def __getattr__(self, attr):
        #allow retrieval of e.g. Ix by calling spin.Ix
        return self.spinOps[attr]


class spinSystem(object):
    def __init__(self, spins, verbose = False, basis = "zeeman", geometry = None):
        """Spins: A list of spins, e.g. ["1H", "1H", "1H", "13C"]"""

        self.spinSystem = []

        self.dimension = 1

        self.geometry = geometry

        self.tol = 1e-15
        
        for s in spins:
            S = Spin(s)
            self.spinSystem.append(S)

            self.dimension = int(self.dimension*(2*S.spin+1))

        maxVal = tuple([s.spin for s in self.spinSystem])
        ranges = [np.arange(-s.spin, s.spin + 1) for s in self.spinSystem]

        #note that we use reverse here so that the matrix representation of singlet
        #spin operators corresponds to the literature.
        elementaryStatesList = reversed(list(itertools.product(*ranges)))
            
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

    def __len__(self):
        return len(self.spinSystem)
    
    def addGeometry(self, predefined = ""):
        self.geometry = spinthon.spin.geometry.Geometry(self, predefined)

        assert(len(self.geometry) == len(self))



        
if __name__ == "__main__":

    methyl = spinSystem(["1H", "1H", "1H"])
    
    methyl13C = spinSystem(["1H", "1H", "1H", "13C"])

    h2o17 = spinSystem(["1H", "1H", "17O"])
    
