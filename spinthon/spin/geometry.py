import numpy as np

class Geometry(object):
    """Represents a molecular geometry. 

    Only the positions of nuclei with non-zero spin are stored.
    
    Public attributes:
    - positions - a list of three-entry numpy arrays, each of which stores coordinates of the respective nucleus
    - name - a descriptive name of the spin system.
    """
    
    def __init__(self, spinSystem, predefined = ""):
        self.positions = []
        self.name = None

        if predefined == "O17water":
            self.name = "O17 Water"
            bondLength = 1.0002e-10
            thetaHOH = 105

            pos1 = np.array([-np.sin(thetaHOH/2), 0 , np.cos(thetaHOH/2)])*bondLength
            pos2 = np.array([np.sin(thetaHOH/2), 0, np.cos(thetaHOH/2)])*bondLength
            pos3 = np.array([0,0,0])*bondLength

            self.positions = [pos1, pos2, pos3]

        elif predefined == "water":
            self.name = "water"
            bondLength = 1.0002e-10
            thetaHOH = 105

            pos1 = np.array([-np.sin(thetaHOH/2), 0 , np.cos(thetaHOH/2)])*bondLength
            pos2 = np.array([np.sin(thetaHOH/2), 0, np.cos(thetaHOH/2)])*bondLength

            self.positions = [pos1, pos2]
            
    def setMiniDNPGeometry(self, e_e_distance = 1e-9, theta = np.pi/4, e_n1_distance = 0.2e-9, theta2 = np.pi/5, phi = np.pi):
            self.name = "miniDNP1"
            
            r1 = e_e_distance
            r2 = e_n1_distance
            
            theta = np.pi/4
            
            pos1 = np.array([0,0,0])
            pos2 = np.array([r1*np.sin(theta), 0, r1*np.cos(theta)])
            pos3 = np.array([r2*np.sin(theta)*np.cos(phi), r2*np.sin(theta)*np.sin(phi), r2*np.cos(theta)])

            self.positions = [pos1, pos2, pos3]
            
    def __getitem__(self, pos):
        return self.positions[pos]


    def __len__(self):
        return len(self.positions)

    def __repr__(self):
        return self.name


    def distance(self, i1, i2):
        return np.linalg.norm(self[i2] - self[i1])
