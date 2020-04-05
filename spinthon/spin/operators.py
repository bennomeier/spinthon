"""Spin operators.

Spin Operators are objects that are defined by a position, i.e. on which spin they act, and a type, i.e. I+, Im, Ix, Iy, Iz.

Spin Operators can be applied to Kets to give new Kets.
It is also possible to get the matrix representation of a spin operator.

Spin Operators can be multiplied.
"""

import numpy as np
from scipy.linalg import expm

from spinthon.basis.vector import Ket


class spinOperator(object):
    def __init__(self, spinSystem, position, which):
        """spin Operator basis class.

        position is an integer that specifies which spin we are talking about.

        which: any of Ix, Iy, Iz, Ip, Im, 1"""

        self.spinSystem = spinSystem
        self.basis = spinSystem.basis

        self.position = position
        self.which = which

        self.matrix = self.getMatrixRepresentation()

    def __mul__(self, other):
        """Presumably this is intended to carry out multiplication.

        It works if the other entity is a Ket."""
        
        if isinstance(other, Ket):
            coeffs = other.coeffs
            ekets = other.ekets

            coeffsNew = np.zeros(np.size(coeffs), dtype = complex)
            retVal = None
            
            if self.which == "Iz":
                for i, c in enumerate(coeffs):
                    coeffsNew[i] = c*ekets[i].Iz(self.position)
                retVal = Ket(coeffsNew, ekets)

            elif self.which == "Ip":
                for i, c in enumerate(coeffs):
                    f, K = ekets[i].Ip(self.position)
                    if K is not None:
                        posNew = ekets.index(K)
                        coeffsNew[posNew] = c*f
                retVal = Ket(coeffsNew, ekets)

            # the following operators are not necessarily correct.
            elif self.which == "Im":
                for i, c in enumerate(coeffs):
                    f, K = ekets[i].Im(self.position)
                    if K is not None:
                        posNew = ekets.index(K)
                        coeffsNew[posNew] += c*f
                retVal = Ket(coeffsNew, ekets)

            elif self.which == "Ix":
                for i, c in enumerate(coeffs):
                    f1, K1 = ekets[i].Ip(self.position)
                    f2, K2 = ekets[i].Im(self.position)

                    if K1 is not None:
                        posNew1 = ekets.index(K1)
                        coeffsNew[posNew1] += c*f1*0.5
                    if K2 is not None:
                        posNew2 = ekets.index(K2)
                        coeffsNew[posNew2] += c*f2*0.5
                retVal = Ket(coeffsNew, ekets)

            elif self.which == "Iy":
                for i, c in enumerate(coeffs):
                    f1, K1 = ekets[i].Ip(self.position)
                    f2, K2 = ekets[i].Im(self.position)

                    if K1 is not None:
                        posNew1 = ekets.index(K1)
                        coeffsNew[posNew1] += complex(0, -c*f1*0.5*1)
                    if K2 is not None:
                        posNew2 = ekets.index(K2)
                        coeffsNew[posNew2] -= complex(0, -c*f2*0.5*1)
                retVal = Ket(coeffsNew, ekets)

            elif self.which == "1":
                retVal = Ket(coeffs, ekets)
                            
            else:
                raise ValueError("Operator Type not supported")
                    
        return retVal


    def getMatrixRepresentation(self):
        """Assume current spin system basis to calculate matrix representation.

        These should compare favorably e.g. with Spin Dynamics, p. 164"""
        
        S = self.spinSystem
        dim = S.dimension
        matRep = np.zeros([dim, dim], dtype = complex)

        for i, b1 in enumerate(S.basis.Bras):
            for j, b2 in enumerate(S.basis.Kets):
                ket = self.__mul__(b2)

                matRep[i,j] += b1*ket

        return matRep


class rotationOperatorSingle(spinOperator):
    def __init__(self, spinSystem, position, beta, phi):
        """a rotation operator for a single spin of the spin system,
        at position position, with flip angle beta and phase phi

        To apply the same rotatition to all spins you need to multiply over all positions.
        """

        Ix = spinOperator(spinSystem, position, "Ix")
        Iy = spinOperator(spinSystem, position, "Iy")

        Ixm = Ix.getMatrixRepresentation()
        Iym = Iy.getMatrixRepresentation()

        tol = spinSystem.tol
        
        matrix =  expm(-1j*beta*(Ixm*np.cos(phi)+Iym*np.sin(phi)))

        matrix.real[abs(matrix.real) < tol] = 0.0
        matrix.imag[abs(matrix.imag) < tol] = 0.0

        self.matrix = matrix
        
class rotationOperatorDirect(spinOperator):
    """ Direct calculation uses a direct product of single spin rotation operator matrix representations,
    as described in Levitt, section 15.8.1. Rotations of a single spin pair"""
    def __init__(self, spinSystem, beta, phi):
        matrix = 1
        for pos in range(len(spinSystem)):
            matrix = np.kron(matrix, spinSystem[pos].rotationOp(beta, phi))

        matrix.real[abs(matrix.real) < spinSystem.tol] = 0.0
        matrix.imag[abs(matrix.imag) < spinSystem.tol] = 0.0
            
        self.matrix = matrix

class rotationOperator(spinOperator):
    """Here, to get the complete rotationOperator we simply use matrix multiplication to multiply the
    rotation operator matrix representations of all spins in the spin system.

    Note that we may want to add a type if we only want to flip say 1H spins in a hetero-nuclear system.

    The approach here should be equivalent to the direct approach, at least in the case of a Zeeman product basis.
    """

    def __init__(self, spinSystem, beta, phi):

        matrix = np.identity(spinSystem.dimension, dtype = complex)
        
        for pos in range(len(spinSystem)):
            matrix = matrix @ rotationOperatorSingle(spinSystem, pos, beta, phi).matrix

        self.matrix = matrix


        
        

    
            
            

        



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



        
