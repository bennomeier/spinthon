"""Spin operators.

This module supplies 5 different kinds of objects

SpinOperatorSingle
These are the elementary operators of which the other Operators are constructed. Examples are I1x, I2z, etc.
Elementary operators may be multiplied with Kets to give new Kets.

SpinOperatorProduct
These are products of elementary operators. Examples are 2*I1x*I2z.
A product is a coefficent, and a list of elementary operators.
Overrides multiplication

SpinOperatorSum
The most general spin operator is a sum of spin operator products. 
Overrides addition and subtraction.

SpinOperator
This is a helper object which creates a SpinOperatorSum and populates it with a single spin operator.

VectorSpinOperator
This is a vector of spinoperators, i.e. S = [Sx, Sy, Sz] with a definition for matrix multiplication.
"""
from array import array
import numpy as np
from scipy.linalg import expm

from spinthon.basis.vector import Ket

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
#SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

def flatten(l):
    return [item for sublist in l for item in sublist]
                                   
class SpinOperatorSingle(object):
    def __init__(self, spinSystem, position, which):
        """spin Operator basis class.

        position is an integer that specifies which spin we are talking about.

        which: any of Ix, Iy, Iz, Ip, Im, 1
        
        """
        

        self.spinSystem = spinSystem
        self.basis = spinSystem.basis

        self.position = position
        self.which = which

        self.matrix = self.getMatrixRepresentation()

    def __mul__(self, other):
        """Presumably this is intended to carry out multiplication.

        It works if the other entity is a Ket."""
        
        if isinstance(other, Ket):
            coeffs = other.coeffs
            ekets = other.ekets

            coeffsNew = np.zeros(np.size(coeffs), dtype = complex)
            retVal = None
            
            if self.which == "Iz":
                for i, c in enumerate(coeffs):
                    coeffsNew[i] = c*ekets[i].Iz(self.position)
                retVal = Ket(coeffsNew, ekets)

            elif self.which == "Ip":
                for i, c in enumerate(coeffs):
                    f, K = ekets[i].Ip(self.position)
                    if K is not None:
                        posNew = ekets.index(K)
                        coeffsNew[posNew] = c*f
                retVal = Ket(coeffsNew, ekets)

            # the following operators are not necessarily correct.
            elif self.which == "Im":
                for i, c in enumerate(coeffs):
                    f, K = ekets[i].Im(self.position)
                    if K is not None:
                        posNew = ekets.index(K)
                        coeffsNew[posNew] += c*f
                retVal = Ket(coeffsNew, ekets)

            elif self.which == "Ix":
                for i, c in enumerate(coeffs):
                    f1, K1 = ekets[i].Ip(self.position)
                    f2, K2 = ekets[i].Im(self.position)

                    if K1 is not None:
                        posNew1 = ekets.index(K1)
                        coeffsNew[posNew1] += c*f1*0.5
                    if K2 is not None:
                        posNew2 = ekets.index(K2)
                        coeffsNew[posNew2] += c*f2*0.5
                retVal = Ket(coeffsNew, ekets)

            elif self.which == "Iy":
                for i, c in enumerate(coeffs):
                    f1, K1 = ekets[i].Ip(self.position)
                    f2, K2 = ekets[i].Im(self.position)

                    if K1 is not None:
                        posNew1 = ekets.index(K1)
                        coeffsNew[posNew1] += complex(0, -c*f1*0.5*1)
                    if K2 is not None:
                        posNew2 = ekets.index(K2)
                        coeffsNew[posNew2] -= complex(0, -c*f2*0.5*1)
                retVal = Ket(coeffsNew, ekets)

            elif self.which == "1":
                retVal = Ket(coeffs, ekets)
                            
            else:
                raise ValueError("Operator Type not supported")
                    
        return retVal
    
    def __repr__(self):
        return   self.which + "{}".format(self.position).translate(SUB)
        #        return    u'I\u208{}O\u209'.format(str(self.position)) + str(self.which)


    def getMatrixRepresentation(self):
        """Assume current spin system basis to calculate matrix representation.

        These should compare favorably e.g. with Spin Dynamics, p. 164"""        
        S = self.spinSystem
        dim = S.dimension
        matRep = np.zeros([dim, dim], dtype = complex)

        for i, b1 in enumerate(S.basis.Bras):
            for j, b2 in enumerate(S.basis.Kets):
                ket = self.__mul__(b2)

                matRep[i,j] += b1*ket

        return matRep


class SpinOperatorProduct(object):
    
    def __init__(self, factor, operators):
        """
        Spin Operator Products are products of elementary spin operators.
        
        factor: a scalar factor
        operators: a list of elementary spin operators

        For example the operator 2 I1x I2y has factor =2 and operators = [I1x, I2y]
        """
        self.factor = factor
        self.operators = operators
        
    def __mul__(self, other):
        """
        Define Multiplication.
        
        Multiplication with a Ket gives a Ket.
        
        Multiplication with a SpinOperatorProduct gives a SpinOperatorProduct.
        
        Multiplication with a number changes the factor.
        """
        # print("Type other: ", type(other))

        
        if isinstance(other, Ket):
            coeffsNew = np.zeros(len(other.coeffs))
            eketsNew = other.ekets
            
            ketNew = Ket(coeffsNew, eketsNew)            
            ketThis = other

            for op in self.operators[::-1]:
                ketThis = op * ketThis    
            
            retVal = self.factor*ketThis
            return retVal
            
        elif isinstance(other, SpinOperatorProduct):
            return SpinOperatorProduct(self.factor*other.factor, self.operators + other.operators)
            
        elif isinstance(other, (int, float)):
            return SpinOperatorProduct(self.factor*other, self.operators)
    
    def __rmul__(self, other):
        return self.__mul__(other)    
            
    def __repr__(self):
        return "{:+}".format(self.factor) + " ".join([o.__repr__() for o in self.operators])
    
    def getMatrixRepresentation(self):
        allMatrices = [o.getMatrixRepresentation() for o in self.operators]
        
        if len(allMatrices) > 1:
            return self.factor*np.linalg.multi_dot(allMatrices)
        else:
            return self.factor*allMatrices[0]
                
        
        
class SpinOperatorSum(object):
    def __init__(self, products):
        """Spin Opeator Sums are sums of spin operator products.
        
        products: a list of spin operator products.
        """
        # sometimes calculations give operators with factor zero. These we want to drop.        
        self.products = [p for p in products if p.factor != 0]
        
    def __mul__(self, other):
        """
        Multiplication proceeds by pairwise multiplication of all products
        """
        if isinstance(other, Ket):            
            coeffsNew = np.zeros(len(other.coeffs))
            eketsNew = other.ekets
            
            ketNew = Ket(coeffsNew, eketsNew)
            for p in self.products:
                ketNew += p*other
            return ketNew
        
        elif isinstance(other, SpinOperatorSum):
            productsNew = []
            for p1 in self.products:
                for p2 in other.products:
                    productsNew.append(p1*p2)
                    
            productsNew = [p1*p2 for p1 in self.products for p2 in other.products]
            return SpinOperatorSum(productsNew)
        
        elif isinstance(other, (SpinOperator)):
            print("It is a Spin operator")
            productsNew = []
            for p1 in self.products:
                for p2 in other.products:
                    productsNew.append(p1*p2)
                    
            productsNew = [p1*p2 for p1 in self.products for p2 in other.products]
            return SpinOperatorSum(productsNew)
        
        elif isinstance(other, (int, float)):
            productsNew = []
            return SpinOperatorSum([SpinOperatorProduct(p.factor*other, p.operators) for p in self.products])
        
        else:
            print(type(other))
            print('Invalid Object for Multiplication with SpinOperatorSum', type(other))
            return NotImplemented
        
    def __rmul__(self, other):
        return self.__mul__(other)
        
    def __add__(self, other):
        #print("Other type: ", type(other))
        assert isinstance(other, (SpinOperatorSum, SpinOperator)), ("In compatible type ", type(other))
        return SpinOperatorSum(self.products + other.products)
    
    def __sub__(self, other):
        assert isinstance(other, (SpinOperatorSum, SpinOperator))
        return self + (-1)*other
    
    def __neg__(self):
        return (-1)*self
    
    def getMatrixRepresentation(self, chop = 0):
        allMatrices = [p.getMatrixRepresentation() for p in self.products]
        s = sum(allMatrices)

        if chop > 0:
            s = s[np.abs(s) < chop]

        return s

    def getEigenValuesAndVectors(self):
        m = self.getMatrixRepresentation()
        eigen_values, eigen_vectors = np.linalg.eigh(m)
        return eigen_values, eigen_vectors

    def getEigenVectors(self):
        m = self.getMatrixRepresentation()
        eigen_values, eigen_vectors = np.linalg.eigh(m)
        return eigen_vectors
        
    
    def __repr__(self):
        return " ".join([p.__repr__() for p in self.products])

class SpinOperator(SpinOperatorSum):
    def __init__(self, spinSystem, position, which):
        """
        SpinOperator - Use this class to generate elementary spin operators.
        
        The class inherits from SpinOperatorSum object, so that operators can be multiplied, and added at will.
        
        Simply inheriting the multiplication routines causes problems with isinstance.
        """
        op = SpinOperatorSingle(spinSystem, position, which)
        self.products = [SpinOperatorProduct(1, [op])]

        
class VectorSpinOperator(object):
    # note the line below is needed STRICTLY,
    # otherwise numpy takes care of multiplication and causes a mess.
    __array_priority__ = 100
    def __init__(self, spinSystem, position):
        self._components = [SpinOperator(spinSystem, position, "Ix"),
                     SpinOperator(spinSystem, position, "Iy"),
                     SpinOperator(spinSystem, position, "Iz")]

    def __iter__(self):
        return iter(self._components)

    def __len__(self):
        return len(self._components)

    def __getitem__(self, index):
        return self._components[index]
        
    def __matmul__(self, other):
        #print("Matmul")
        #print("Self: ", self)
        print("Other: ", other.__class__.__mro__)
        #0print("Type Other: ", type(other).__name__)
        
        if isinstance(other, np.ndarray):
            # Note that here we just assume that other is a matrix."""
            # print("Shape: ", other.shape)
            if other.shape == (3,3):
                # print("Other is 3x3 Matrix.")
                c1 = other[:,0]@self._components
                # print("c1: ", c1)
                c2 = other[:,1]@self._components
                c3 = other[:,2]@self._components
                return VectorSpinOperatorGeneral([c1, c2, c3])
            elif other.shape ==(3,):
                # print("Multiplying with array.")
                return other[0]*self[0] + other[1]*self[1] + other[2]*self[2]
        
        # note that below we would normally use isinstance(), however 
        # for some arcane reason isinstance(other, VectorSpinOperator) does not
        # seem to evaluate correctly.
        
        elif type(other).__name__ == "VectorSpinOperator" or type(other).__name__ == "VectorSpinOperatorGeneral":
            #print("Now we should return a simple spinOperator")
            res = self[0]*other[0] + self[1]*other[1] + self[2]*other[2]
            # print(res)
            return res
        
        else:
            print("Something went wrong.")
            print("other: ", other)
            print("self: ", self)
            print(other.__class__.__mro__)
            print("type (other2): ", type(other).__name__)
            print(type(other) ==  "<class 'spinthon.spin.operators.VectorSpinOperator'>")
            print("type (self): ", type(self))
            return NotImplemented
    

                
    #def __array_wrap__(self, result):
    #    """This is required so that multiplication is not carried out by numpy."""
     #   return VectorSpinOperatorGeneral(self._components)  # can add other attributes of self as constructor

    def __rmatmul__(self, other):
        """ Right Multiplication. Self appears on the right,
        other on the left. Self is a spin vector, other can be an np.array"""
        if isinstance(other, np.ndarray):
            if other.shape == (3,3):
                # print("Other is 3x3 Matrix.")
                c1 = other[:,0]@self._components
                # print("c1: ", c1)
                c2 = other[:,1]@self._components
                c3 = other[:,2]@self._components
                return VectorSpinOperatorGeneral([c1, c2, c3])
                        
            elif other.shape ==(3,):
                # print("Multiplying with array.")
                opR = other[0]*self[0] + other[1]*self[1] + other[2]*self[2]
                return opR
                
            
            else:
                print("Rmatmul", other)
                print("self: ", self)
                print("other: ", other)
                print("type (other2): ", type(other))

                print("THIS IS NOT IMPLEMENTED.")
        
                return NotImplemented

    def __repr__(self):
        return "[" + ", ".join([p.__repr__() for p in self._components]) + "]"

    

class VectorSpinOperatorGeneral(VectorSpinOperator):
    def __init__(self, components):
        self._components = components