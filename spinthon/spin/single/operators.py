import numpy as np
from scipy.linalg import expm
#from scipy.constants import k, hbar

def getSpinOperators(I, retOnly=""):
    """getSpinOperators(I) returns a dictionary containing entries Ix, Iy, Iz, Ip, Im
    for any integer or half inter spin I"""

    ops = {}
    dim = int(2*I + 1)
    zComponent = np.arange(I, -I-1, -1)
    Iz = np.diag(zComponent)

    Iplus = np.zeros([dim,dim])
    Iminus = np.zeros([dim,dim])

    for k in range(int(dim-1)):
        m = zComponent[k+1]
        Iplus[k,k+1] = np.sqrt(I*(I+1) - m*(m+1))
        m = zComponent[k]
        Iminus[k+1,k] = np.sqrt(I*(I+1) -m*(m-1))

    Ix = 0.5*(Iplus + Iminus)
    Iy = -0.5j*(Iplus - Iminus)

    ops = {"Ix": Ix, "Iy": Iy, "Iz" : Iz, "Iplus": Iplus, "Iminus": Iminus}
    if retOnly == "":
        retVal = ops
    else:
        retVal = ops[retOnly]
    return retVal
        

def getSingleSpinRotationOperator(I, phi, beta):
    spinOps = getSpinOperators(I)
    Ix = spinOps["Ix"]
    Iy = spinOps["Iy"]
    return expm(-1j*beta*(Ix*np.cos(phi) + Iy*np.sin(phi)))


def getRotationOperator(matrix):
    return expm(-1j*matrix)

def getRotationOperatorList(list):
    return [expm(-1j*element) for element in list]

#def getEqDensityMatrix(I,T,nu):
#    if I != 0.5:
#        print "getDensityMatrix currently worksOnly for spin 1/2"
#        
#    B = hbar*2*np.pi*nu/(k*T)
#    print "The Boltzmann factor is {0:1e}".format(B)
 #   rhoEq = np.matrix([[0.5 + 0.25*B, 0],[0, 0.5 - 0.25*B]])
#    if I == 0.5:
#        return rhoEq


if __name__ == "__main__":
    print("spinOperators Test")
    ops = getSpinOperators(1.5)

    for i in ops:
        print()
        print(i)
        print(ops[i])
