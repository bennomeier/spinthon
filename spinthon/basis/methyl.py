import itertools
import numpy as np
#import copy

from zeeman import zeemanProductBasis

class methyl_A_E(zeemanProductBasis):
    def __init__(self, spinSystem):
        zeemanProductBasis.__init__(self, spinSystem)

        #we have to refactor the first three spins
        
        #assume the singlet-triplet quantum numbers come first.
        #hence part all the kets into four sets of the same length

        #the first ones are b,b,b,x, followed by b,b,a,x, ..., a, a, a, x 

        self.name = "Methyl A and E states basis"
        
        dim = spinSystem.dimension

        # this is the multiplicity of each of the a and e states.
        # if there is another spin 1/2 in the system, the multiplicity of all states is 2.
        multiplicity = dim/8 
        

        self.basisNew = []

        #A states, starting with b,b,b, total spin -3/2 
        self.basisNew.extend(self.basis[:multiplicity])

        #next A states with total spin -1/2
        for k in range(multiplicity):
            self.basisNew.append(1/np.sqrt(3)*(self.basis[4*multiplicity + k] +
                                 self.basis[1*multiplicity + k] +
                                 self.basis[2*multiplicity + k]))

        #next A states with total spin 1/2
        for k in range(multiplicity):
            self.basisNew.append(1/np.sqrt(3)*
                                 (self.basis[6*multiplicity + k] +
                                 self.basis[5*multiplicity + k] +
                                 self.basis[3*multiplicity + k]))

        self.basisNew.extend(self.basis[7:7+multiplicity])

        #E states, start with Ea, -1/2
        epsilon = np.exp(1j*2*np.pi/3)
        for k in range(multiplicity):
            self.basisNew.append(1/np.sqrt(3)*
                                 (self.basis[4*multiplicity + k] +
                                  epsilon*self.basis[1*multiplicity + k] +
                                  np.conjugate(epsilon)*self.basis[2*multiplicity + k]))

        #Ea, total spin 1/2
        for k in range(multiplicity):
            self.basisNew.append(1/np.sqrt(3)*
                                 (self.basis[6*multiplicity + k] +
                                  epsilon*self.basis[5*multiplicity + k] +
                                  np.conjugate(epsilon)*self.basis[3*multiplicity + k]))

        #Eb states are like Ea with conjugation swapped; -1/2
        for k in range(multiplicity):
            self.basisNew.append(1/np.sqrt(3)*
                                 (self.basis[4*multiplicity + k] +
                                  np.conjugate(epsilon)*self.basis[1*multiplicity + k] +
                                  epsilon*self.basis[2*multiplicity + k]))

        #Eb, total spin 1/2
        for k in range(multiplicity):
            self.basisNew.append(1/np.sqrt(3)*
                                 (self.basis[6*multiplicity + k] +
                                  np.conjugate(epsilon)*self.basis[5*multiplicity + k] +
                                  epsilon*self.basis[3*multiplicity + k]))
            
            
        self.basis = self.basisNew

