import numpy as np

def Rx(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s, c]])

def Ry(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, 0, s],
                    [0, 1, 0],
                    [-s, 0, c]])


def Rz(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, -s, 0],
                    [s, c, 0],
                    [0, 0, 1]])

def R(alpha, beta, gamma, tol = 1e-16):
    """Euler Rotation Matrix using the zyz convention,
    cf. Malcolm Levitt, Spin Dynamics, 2nd Edition, p. 601"""
    
    ca, cb, cg = np.cos(alpha), np.cos(beta), np.cos(gamma)
    sa, sb, sg = np.sin(alpha), np.sin(beta), np.sin(gamma)

    m = np.array([[ca*cb*cg - sa*sg, -sa*cg - ca*cb*sg, ca*sb],
                  [sa*cb*cg + ca*sg, ca*cg - sa*cb*sg, sa*sb],
                  [-sb*cg, sb*sg, cb]])

    #m[np.abs(m) < tol] = 0
    return m


def tensorTransform(alpha, beta, gamma, t):
    # using the notation by Levitt (Spin Dynamics),
    # R as defined above transforms an A system vector into a B system vector
    # 

    # for passive transforms from axis system A (e.g. chemical shift principal axis frame) to
    # B (e.g. lab frame), we have

    
    Rleft = R(alpha, beta, gamma)
    # for rotation matrices, the transpose is also the invere.
    Rright =  np.transpose(Rleft)

    return Rleft@t@Rright


if __name__ == "__main__":
    alpha = 0.2
    beta = 0.9
    gamma = 1.3

    print("Euler Matrix: ", R(alpha, beta, gamma))
    print("Rz(gamma)*Ry(beta)*Rz(alpha):", Rz(gamma)@Ry(beta)@Rz(alpha))

    print("Inverse Check: ", R(alpha, beta, gamma)@R(-gamma, -beta, -alpha))

    print("Check inverse is transpoes: ", R(alpha, beta, gamma)@np.transpose(R(alpha, beta, gamma)))

                    
