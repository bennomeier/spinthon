import numpy as np

def d(m,n,beta):
    """Return the reduced Wigner rotation matrices
    d^2_{mn}(beta)

    m, n are in [-2,...,2]
    beta: Euler angle, in radians
    """

    cb = np.cos(beta)
    sb = np.sin(beta)
    
    # the reduced wigner rotation matrix elements, cf. e.g. Protein NMR Spectroscopy, Second Edition, p.105
    matrix = np.array([[np.cos(beta/2.)**4, 1./2*sb*(cb + 1), np.sqrt(3./8)*sb**2, -1./2*sb*(cb -1), np.sin(beta/2.)**4],
                       [-1./2*sb*(cb + 1), 1./2*(2*cb -1)*(cb + 1), np.sqrt(3./2)*sb*cb, -1./2*(2*cb -1)*(cb - 1), -1./2*sb*(cb - 1)],
                       [np.sqrt(3./8)*sb**2, - np.sqrt(3/2)*sb*cb, 1./2*(3*cb**2 - 1), np.sqrt(3./2)*sb*cb, np.sqrt(3./8)*sb**2],
                       [1./2*sb*(cb -1), - 1./2*(2*cb - 1)*(cb - 1), - np.sqrt(3./2)*sb*cb, 1./2*(2*cb - 1)*(cb + 1), 1./2*sb*(cb + 1)],
                       [np.sin(beta/2.)**4, 1./2*sb*(cb-1), np.sqrt(3./8)*sb**2, -1./2*sb*(cb + 1), np.cos(beta/2.)**4]])

    print(matrix)

    return matrix[m+2, n+2]

    
def D(m,n,alpha,beta,gamma):
    """Returns the Wigner Rotation Matrix Element D_{mn}

    Arguments:
    m: row index
    n: column index
    alpha: Euler Angle, in radians
    beta: Euler Angle, in radians
    gamma: Euler Angle, in radians
    """
    
    return d(m,n,beta)*np.exp(-1j*m*alpha - 1j*n*gamma)

if __name__ == "__main__":
    from euler import tensorTransform

    tensor = np.diag([0.5, 1.3, 2.2])
    print(tensorTransform(1, 1.3, 2, tensor))
    
    tensorPrime = np.zeros(np.shape(tensor))
    
    
