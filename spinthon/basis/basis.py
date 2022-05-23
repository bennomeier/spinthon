import numpy as np
import itertools

from spinthon.basis.vector import Ket, Bra
from spinthon.spin.operators import SpinOperator

class Basis(object):
    def __init__(self, spinSystem, kets, name = "User_Defined_Basis"):
        """A Basis is just a list of kets and corresponding bras."""
        self.spinSystem = spinSystem
        self.name = name
        self.Kets = kets
        self.Bras = [k.getBra() for k in self.Kets]

    def __getitem__(self, pos):
        """support indexing of the basis."""
        return self.Kets[pos]

    def __repr__(self):
        ans = [k.__repr__() for k in self.Kets]
        return "\n".join(ans)

    def getAllowedTransitions(self, position, change):
        return None

    def rebase(self, eigenvectors):
        """ to calculate the new basis we rotate each 
        basis vector with the matrix of eigenvectors."""
        
        M = eigenvectors

        #M = [ev/np.sqrt(np.linalg.norm(ev)) for ev in eigenvectors]
        
        ekets = self.Kets[0].ekets
        KetsNew = []
        for i in range(len(ekets)):            
            thisKet = Ket(M@self.Kets[i].coeffs, ekets)
            thisKet = 1/np.linalg.norm(thisKet.coeffs)*thisKet

            KetsNew.append(thisKet)
        # self.Kets = [Ket(M[i], ekets) for i in range(len(ekets))]
        self.Kets = KetsNew
        self.Bras = [k.getBra() for k in self.Kets]
        

    def checkOrthoNormal(self):
        Op = SpinOperator(self.spinSystem, 0, "1")
        print(Op.getMatrixRepresentation(chop = 0))

    def getAllowedTransitions(self, changes):
        """Calculate all the allowed transitions. 
        
        changes:   the allowed changes. For example, for a two spin system, if we are interested in changes in the first spin, these would be
        
        [[1, 0], [-1, 0]]
        """
        combinations = itertools.combinations(self.Kets, 2)
        print(combinations)