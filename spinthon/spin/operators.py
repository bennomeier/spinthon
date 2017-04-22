"""Calculate spin operators for a given basis"""

import numpy as np



# a spin operator should be an object. It can only have a matrix representation with respect to
# a certain base. The base should hence be stored with the matrix representation.

# operator representations with the same basis can be added, multiplied, etc.


class spinOperator(object):
    def __init__(self, spinSystem, position, which):
        """spin Operator basis class.

        which: any of Ix, Iy, Iz, Ip, Im"""

        self.spinSystem = spinSystem
        self.basis = spinSystem.basis

        pos = self.position

        



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



        
