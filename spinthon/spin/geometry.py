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

            
    def __getitem__(self, pos):
        return self.positions[pos]


    def __len__(self):
        return len(self.positions)

    def __repr__(self):
        return self.name


    def distance(self, i1, i2):
        return np.linalg.norm(self[i2] - self[i1])
