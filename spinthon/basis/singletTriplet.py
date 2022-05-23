"""This module provides a class that defines a basis
suitable for singlet/triplet spin states."""

import numpy as np


from spinthon.basis.zeeman import ZeemanProductBasis
from spinthon.spin.system import SpinSystem


class SingletTripletBasis(ZeemanProductBasis):
    """singletTripletBasis provides a new basis for a given spinSystem.
    Its init routine calls the parent class constructor, zeemanProductBasis.
    The singlet/triplet pair comprises the first two spins in the spin system.
    New linear combinations of the elementary product kets are formed.

    The basis can be divided into four subspaces.

    The first basis vectors up to dim/4 correspond to the singlet state
    of the first two spins in the sytem.
    The second batch (i.e. from dim/4:2*dim/4) corresponds to the T_{+1} states
    , then follow the T_{0} states and ultimately the T_{-1} states.
    """
    def __init__(self, spinSystem):
        ZeemanProductBasis.__init__(self, spinSystem)

        # by definition the singlet-triplet quantum numbers come first.
        # hence part all the kets into four sets of the same length

        # the first ones are +1/2, +1/2,x
        # the second ones are 1/2, -1/2, x
        # the third ones are -1/2, 1/2, x
        # the fourth ones are -1/2, -1/2, x

        self.name = "Singlet Triplet Basis"
        self.spinSystem = spinSystem

        # hence first and fourths onces are triplets,
        # second and third ones need to be mixed

        dim = spinSystem.dimension

        start = int(dim/4)
        stop = int(2*dim/4)

        for k in range(start, stop):
            bFix = self.Kets[k]
            self.Kets[k] = 1/np.sqrt(2)*(self.Kets[k + start] + self.Kets[k])
            self.Kets[k + start] = 1/np.sqrt(2)*(self.Kets[k+start] - bFix)

        # now lets rearange:

        self.Kets = (self.Kets[2*start:3*start] +
                     self.Kets[0:2*start] + self.Kets[3*start:])
        self.Bras = [k.getBra() for k in self.Kets]

        spinSystem.basis = self


if __name__ == "__main__":
    print("Singlet Triplet Basis Test")
    S = SpinSystem(["1H", "1H", "13C"])
    spinSystem.basis = SingletTripletBasis(S)
    print("Basis: " + spinSystem.basis.name)

    for k in S.basis:
        print(k)
