import numpy as np

class BasisKet(object):
    """Basis Kets are immutable objects. 
    They correspond to the zeeman product basis and are stored as tuples.

    This class also supplies methods to apply a basic spin operator to the i-th spin."""
    
    def __init__(self, value, maxVal):
        self.maxVal = maxVal
        self.value = value


    def Iz(self, pos):
        M = self.value[pos]
        return m

    def Ip(self, pos):
        M = self.value[pos]
        I = self.maxVal[pos]

        f = 0
        returnKet = None
        
        if M < I:
            newValue = self.value
            newValue[pos] = M+1
            
            f = np.sqrt(I*(I+1) - M*(M+1))
            returnKet = BasisKet(newValue, self.maxVal)

        return f, returnKet

    def Im(self, pos):
        M = self.value[pos]
        I = self.maxVal[pos]

        f = 0
        returnKet = None
        
        if M > -I:
            newValue = self.value
            newValue[pos] = M-1

            f = np.sqrt(I*(I+1) - M*(M-1))
            returnKet = BasisKet(newValue, self.maxVal)
        
        return f, returnKet
            

    def __repr__(self):
        ketString = "|"
        for t in self.value:
            if len(ketString) > 1:
                ketString += ","
            ketString += str(t)
        ketString += ">"

        return ketString
        

class Ket(object):
    def __init__(self, coeffs, basiskets):
        """A Ket is simply a linear combination of basiskets of the spinsystem.

        Note that in spinthon for any ket the full list of coefficients is stored.
        The length of this list corresponds to the spin system dimension."""
        
        self.coeffs = coeffs
        self.basiskets = basiskets


    def __add__(self, other):
        assert self.basiskets == other.basiskets
        
        return Ket(self.coeffs + other.coeffs, self.basiskets)

    def __sub__(self, other):
        assert self.basiskets == other.basiskets
        
        return Ket(self.coeffs - other.coeffs, self.basiskets)

    def __rmul__(self, factor):
        return Ket(self.coeffs*factor, self.basiskets)

    
    def __repr__(self):
        output = ""
        for k in range(len(self.basiskets)):
            if self.coeffs[k] != 0:


                if len(output) > 0:
                    output += " + "
                
                output += "{:.2f} x".format(self.coeffs[k]) + str(self.basiskets[k]) #+ str(self.kets[k])

        return output

    def Iz(self, pos):
        return np.sum([self.coeffs[k]*self.basiskets[k].Iz(pos) for k in range(len(self.basiskets))])


    def Ip(self, pos):
        """Use the basisket Ip to shift all the basiskets. Now the trick is to know where the coefficients need to go."""

        for k in range(len(self.basiskets)):
            if self.coeffs[k] != 0:
                pass
                
        
