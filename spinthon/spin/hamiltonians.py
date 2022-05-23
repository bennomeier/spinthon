"""Hamiltonians

Hamiltonians are really ony Spin Operator Sums. 

This module simplifies the construction of Hamiltonians by providing standard Hamilltonians.

The Hamiltonian class also allows the extraction of energies for a given transaction.
"""

import numpy as np

from spinthon.basis.vector import Ket
from spinthon.spin.operators import SpinOperatorSum, SpinOperator

class Hamiltonian(SpinOperatorSum):
    def __init__(self):
        self.products = []
    
    
    def diagonalize(self, system, rebase = True):
        Hmat = self.getMatrixRepresentation()
        energies, eigenVectors = np.linalg.eig(Hmat)
        
        if rebase:
            for i in range(system.dim):
                system.basis.Kets[i].coeffs = eigenVectors[i]
                system.basis.Bras[i].coeffs = eigenVectors[i]
        
        return np.real(energies)
        
    
    
class ElectronZeeman(Hamiltonian):
    def __init__(self, system, position, B, g):
        """Anisotropic Zeeman Interaction Zeeman Interaction
        
        system: spin system
        p: position of the involved spin
        B: magnetic field (3 component numpy array)
        g: g tensor (3x3 numpy array)
        """
        
        Sx = SpinOperator(system, position, "Ix")
        Sy = SpinOperator(system, position, "Iy")
        Sz = SpinOperator(system, position, "Iz")
        
        H = B@g@[Sx, Sy, Sz]
        self.products = H.products
        
class Hyperfine(Hamiltonian):
    def __init__(self, system, posE, posN, A):
        """Anisotropic Hyperfine Interaction
        
        system: spin system
        posE: position of the electron
        posN: position of the nucleus
        A: hyperfine couping tensor (3x3 numpy array)"""
        
        Sx = SpinOperator(system, posE, "Ix")
        Sy = SpinOperator(system, posE, "Iy")
        Sz = SpinOperator(system, posE, "Iz")
        
        Ix = SpinOperator(system, posN, "Ix")
        Iy = SpinOperator(system, posN, "Iy")
        Iz = SpinOperator(system, posN, "Iz")
        
        
        H = [Sx, Sy, Sz]@A@[Ix, Iy, Iz]
        self.products = H.products
        
        
