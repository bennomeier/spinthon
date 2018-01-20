"""Vector Module

The vector module defines elementary bras and kets and their behaviour upon multiplication, inner products, etc.

The module also defines Bras and Kets which are linear combinations of elementary bras and kets. Bras and Kets also comprise the basis of any spin system."""

import numpy as np

class ElementaryKet(object):
    """Elementary Kets are immutable objects. 
    They correspond to the Zeeman product basis kets and are stored as tuples.

    This class also supplies methods to apply a basic spin operator to the i-th spin."""
    
    def __init__(self, value, maxVal, repr = "plusMinus"):
        self.maxVal = maxVal
        self.value = value
        self.repr = repr
        if self.repr == "plusMinus":
            self.reprDict = {0.5: "+", -0.5: "-"}
        elif self.repr == "alphaBeta":
            self.reprDict = {0.5: "α", -0.5: "β"}
        
    def __eq__(self, other):
        return self.value == other.value

    def getBra(self):
        return ElementaryBra(self.value, self.maxVal)

    def Iz(self, pos):
        M = self.value[pos]
        return M

    def Ip(self, pos):
        M = self.value[pos]
        I = self.maxVal[pos]

        f = 0
        returnKet = None
        
        if M < I:
            newValue = list(self.value)
            newValue[pos] = M+1
            
            f = np.sqrt(I*(I+1) - M*(M+1))
            returnKet = ElementaryKet(tuple(newValue), self.maxVal)

        return f, returnKet

    def Im(self, pos):
        M = self.value[pos]
        I = self.maxVal[pos]

        f = 0
        returnKet = None
        
        if M > -I:
            newValue = list(self.value)
            newValue[pos] = M-1

            f = np.sqrt(I*(I+1) - M*(M-1))
            returnKet = ElementaryKet(tuple(newValue), self.maxVal)
        
        return f, returnKet
            

    def __repr__(self):
        ketString = "|"
        for i, t in enumerate(self.value):
            if self.repr != "numeric":
                ketString += self.reprDict[t]
            else:
                if len(ketString) > 1:
                    ketString += ","
                ketString += str(t)
        ketString += ">"

        return ketString


class ElementaryBra(object):
    """Elementary Bras are immutable objects, corresponding to Zeeman product states.
   
    Multiplication with an elementary ket (inner product) yields 0 or 1.
    """
    def __init__(self, value, maxVal, repr = "alphaBeta"):
        self.maxVal = maxVal
        self.value = value
        self.repr = repr
        self.alphaBetaDict = {0.5: "+", -0.5: "-"}

    def __mul__(self, other):
        """Multiplication with a Ket from the right is the inner product."""
        assert(isinstance(other, ElementaryKet))
        retValue = 0
        if other.value == self.value:
            retValue = 1
        return retValue

    def __repr__(self):
        ketString = "<"
        for i, t in enumerate(self.value):
            if self.repr == "alphaBeta":
                ketString += self.alphaBetaDict[t]
            else:
                if len(ketString) > 1:
                    ketString += ","
                ketString += str(t)
        ketString += "|"

        return ketString

    
class Ket(object):
    def __init__(self, coeffs, ekets):
        """A Ket is simply a linear combination of elementary kets of the spinsystem.

        Note that in spinthon for any ket the full list of coefficients is stored.
        The length of this list corresponds to the spin system dimension.

        The ekets are the elementary zeeman product kets, irrespective of the chosen base.
        """
         
        self.coeffs = coeffs
        self.ekets = ekets


    def __add__(self, other):
        assert self.ekets == other.ekets
        
        return Ket(self.coeffs + other.coeffs, self.ekets)

    def __sub__(self, other):
        assert self.ekets == other.ekets
        
        return Ket(self.coeffs - other.coeffs, self.ekets)

    def __rmul__(self, factor):
        return Ket(self.coeffs*factor, self.ekets)

    
    def __repr__(self):
        output = ""
        for k in range(len(self.ekets)):
            if self.coeffs[k] != 0:


                if len(output) > 0:
                    output += " + "
                
                output += "{:.2f} x".format(self.coeffs[k]) + str(self.ekets[k])

        return output

    def getBra(self):
        return Bra(self.coeffs, [b.getBra() for b in self.ekets])

    def Iz(self, pos):
        return np.sum([self.coeffs[k]*self.ekets[k].Iz(pos) for k in range(len(self.ekets))])


    def Ip(self, pos):
        """Use the eket Ip to shift all the ekets. Now the trick is to know where the coefficients need to go."""

        for k in range(len(self.ekets)):
            if self.coeffs[k] != 0:
                pass


class Bra():
    """A bra is a linear combination of elementary bras."""
    
    def __init__(self, coeffs, ebras):
        """A Bra is a linear combination of ebras of the spinsystem.

        Note that in spinthon for any ket the full list of coefficients is stored.
        The length of this list corresponds to the spin system dimension."""
        
        self.coeffs = coeffs
        self.ebras = ebras


    def __add__(self, other):
        assert self.ebras == other.ebras
        
        return Bra(self.coeffs + other.coeffs, self.ebras)

    def __sub__(self, other):
        assert self.ebras == other.ebras
        
        return Bra(self.coeffs - other.coeffs, self.ebras)

    def __mul__(self, other):
        """Multiplication with a Ket from the right is the inner product."""
        assert(isinstance(other, Ket))        

        sum = 0
        for k in range(len(self.coeffs)):
            sum += np.conj(self.coeffs[k])*other.coeffs[k]
        
        return sum
    
    
    def __rmul__(self, factor):
        return Bra(self.coeffs*factor, self.ebras)

    
    def __repr__(self):
        output = ""
        for k in range(len(self.ebras)):
            if self.coeffs[k] != 0:


                if len(output) > 0:
                    output += " + "
                
                output += "{:.2f} x".format(self.coeffs[k]) + str(self.ebras[k])

        return output   
        

        
