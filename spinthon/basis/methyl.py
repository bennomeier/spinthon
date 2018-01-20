import numpy as np

from spinthon.basis.zeeman import zeemanProductBasis
from spinthon.spin.system import spinSystem


class methyl_A_E(zeemanProductBasis):
    def __init__(self, spinSystem):
        zeemanProductBasis.__init__(self, spinSystem)

        # we have to refactor the first three spins
        # assume the singlet-triplet quantum numbers come first.
        # hence part all the kets into four sets of the same length
        # the first ones are b,b,b,x, followed by b,b,a,x, ..., a, a, a, x

        self.name = "Methyl A and E states basis"
        dim = spinSystem.dimension

        # this is the multiplicity of each of the a and e states.
        # if there is another spin 1/2 in the system, the
        # multiplicity of all states is 2.
        multiplicity = int(dim/8)
        self.KetsNew = []

        # A states, starting with a,a,a, total spin 3/2
        self.KetsNew.extend(self.Kets[:multiplicity])

        # next A states with total spin 1/2
        for k in range(multiplicity):
            self.KetsNew.append(1/np.sqrt(3)*(self.Kets[4*multiplicity + k] +
                                              self.Kets[1*multiplicity + k] +
                                              self.Kets[2*multiplicity + k]))

        # next A states with total spin - 1/2
        for k in range(multiplicity):
            self.KetsNew.append(1/np.sqrt(3) *
                                (self.Kets[6*multiplicity + k] +
                                 self.Kets[5*multiplicity + k] +
                                 self.Kets[3*multiplicity + k]))

        self.KetsNew.extend(self.Kets[7*multiplicity:8*multiplicity])

        # E states, start with Ea, -1/2
        epsilon = np.exp(1j*2*np.pi/3)
        for k in range(multiplicity):
            self.KetsNew.append(1/np.sqrt(3) *
                                (self.Kets[4*multiplicity + k] +
                                 epsilon*self.Kets[1*multiplicity + k] +
                                 np.conjugate(epsilon) *
                                 self.Kets[2*multiplicity + k]))

        # Ea, total spin 1/2
        for k in range(multiplicity):
            self.KetsNew.append(1/np.sqrt(3) *
                                (self.Kets[6*multiplicity + k] +
                                 epsilon*self.Kets[5*multiplicity + k] +
                                 np.conjugate(epsilon) *
                                 self.Kets[3*multiplicity + k]))

        # Eb states are like Ea with conjugation swapped; -1/2
        for k in range(multiplicity):
            self.KetsNew.append(1/np.sqrt(3) *
                                (self.Kets[4*multiplicity + k] +
                                 np.conjugate(epsilon) *
                                 self.Kets[1*multiplicity + k] +
                                 epsilon*self.Kets[2*multiplicity + k]))

        # Eb, total spin 1/2
        for k in range(multiplicity):
            self.KetsNew.append(1/np.sqrt(3) *
                                (self.Kets[6*multiplicity + k] +
                                 np.conjugate(epsilon) *
                                 self.Kets[5*multiplicity + k] +
                                 epsilon*self.Kets[3*multiplicity + k]))
        self.Kets = self.KetsNew

        spinSystem.basis = self


if __name__ == "__main__":
    print("Methyl Basis Test")
    S = spinSystem(["1H", "1H", "1H", "13C"])
    spinSystem.basis = methyl_A_E(S)
    print("Basis: " + spinSystem.basis.name)

    for k in S.basis:
        print(k)
