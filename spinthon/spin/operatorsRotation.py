
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



        
