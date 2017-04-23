"""Spin operators.

Spin Operators are objects that are defined by a position, i.e. on which spin they act, and a type, i.e. I+, Im, Ix, Iy, Iz.

Spin Operators can be applied to Kets to give new Kets.
It is also possible to get the matrix representation of a spin operator.

Spin Operators can be multiplied.
"""

import numpy as np

from spinthon.basis.vector import Ket

class spinOperator(object):
    def __init__(self, spinSystem, position, which):
        """spin Operator basis class.

        which: any of Ix, Iy, Iz, Ip, Im"""

        self.spinSystem = spinSystem
        self.basis = spinSystem.basis

        self.position = position
        self.which = which

    def __mul__(self, other):
        if isinstance(other, Ket):
            coeffs = other.coeffs
            ekets = other.ekets

            coeffsNew = np.zeros(np.size(coeffs))
            retVal = None
            
            if self.which == "Iz":
                for i, c in enumerate(coeffs):
                    coeffsNew[i] = c*ekets[i].Iz(self.position)
                retVal = Ket(coeffsNew, ekets)

            if self.which == "Ip":
                for i, c in enumerate(coeffs):
                    f, K = ekets[i].Ip(self.position)
                    if K is not None:
                        posNew = ekets.index(K)
                        coeffsNew[posNew] = c*f
                retVal = Ket(coeffsNew, ekets)
                    
        return retVal


    def getMatrixRepresentation(self):
        """Assume current spin system basis to calculate matrix representation.

        These should compare favorably e.g. with Spin Dynamics, p. 164"""
        
        S = self.spinSystem
        dim = S.dimension
        matRep = np.zeros([dim, dim])

        for i, b1 in enumerate(S.basis.Bras):
            for j, b2 in enumerate(S.basis.Kets):
                ket = self.__mul__(b2)

                matRep[i,j] += b1*ket

        return matRep
                
            
            

        



def getSpinOperators(basis, position, only = ""):
    """return spin operators for the spin at position. 

    with only = "", a dictionary with entries Ix, Iy, Iz, Ip, Im is returned.
    
    if only is a key of that dictionary, only the matrix representation of the specified operator is returned."""


    Iz = np.diag([b.Iz(position) for b in basis])

    #Iplus and Iminus are defined as
    #I+ |I,M> = [I(I+1) - M(M+1)]^0.5 |I, M+1>
    #I- |I,M> = [I(I+1) - M(M-1)]^0.5 |I, M-1>
    
    #Iplus

    return Iz



        
