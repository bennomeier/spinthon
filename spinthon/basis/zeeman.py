#this is required for python > 3.4
#from importlib import reload

import numpy as np

from spinthon.basis.vector import Ket, Bra
from spinthon.basis.basis import Basis

class ZeemanProductBasis(Basis):
    def __init__(self, spinSystem):
        """Elementary zeeman product basis.

        Parent class to all other bases.
        
        Any basis has fields
        
        b.Kets
        b.Bras

        These are linear combinations of the elementary kets / bras of the spin system.

        The __getitem__ method will return the Ket.
        """
        self.spinSystem = spinSystem
        self.name = "Zeeman Product Basis"

        S = spinSystem

        # for the zeeman basis the coefficients are just one in one position,
        # and zeros elsewhere        
        self.Kets = []


        for k in range(S.dimension):
            coefficients = np.zeros(S.dimension)
            coefficients[k] = 1

            #note that it is redundant to store the kets
            self.Kets.append(Ket(coefficients, S.ekets))

        self.Bras = [k.getBra() for k in self.Kets]
        
