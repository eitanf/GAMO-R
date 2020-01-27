"""
Author: Hrishee Shastri
May 2019

A Test function interface 

"""
import math
import random

class TestFn:
    """
    encapsulates all data needed to interact with a given test function for GA optimization

    name -- a string description of the function (e.g. "Parabola with noise")
    formula -- a real valued function that is scalar valued in its output and vector valued in its input   
    dim -- dimension of the input space R^dim
    """
    def __init__(self, name, formula, dimension):
        self._name = name
        self._f = formula 
        self._n = dimension

    def eval(self, vector):
        """
        evaluates the function with a given real valued vector. The vector can be a tuple or a list
        """
        if len(vector) != self._n:
            raise ValueError("Input dimensions don't match")
        return self._f(vector)

    def get_input_dimension(self):
        return self._n

    def __str__(self):
        return self._name




############ Benchmarks: some mathematical functions for optimization ###############


# Rosenbrock's saddle function in 2 dimensions
f2 = TestFn("Rosenbrock's Saddle", lambda X: ((1-X[0])**2)+100*((X[0]**2)-X[1])**2, dimension=2)

# Beale's function in 2 dimensions
BEALEf = TestFn("Beale function", lambda X: ((1.5-X[0]+X[0]*X[1])**2) + ((2.25-X[0]+X[0]*X[1]**2)**2) + ((2.625-X[0]+X[0]*X[1]**3)**2), dimension=2)

# Parabola in 3 dimensions
f1 = TestFn("Parabola", lambda X: sum([x_i**2 for x_i in X]), dimension=3)

# Step function in 5 dimensions
f3 = TestFn("Step function", lambda X: sum([math.floor(X[i]) for i in range(len(X))]), dimension=5)

# Quartic with noise in 30 dimensions
f4 = TestFn("Quartic with noise", lambda X: sum([i*(X[i]**4) for i in range(len(X))]) + random.gauss(mu=0,sigma=1), dimension=30)

# Shekel's foxholes in 2 dimension
def shekel(X):
    A_1 = [-32, -16, 0, 16, 32]*5
    A_2 = [-32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 0, 0, 0, 0, 0, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32,]
    A = [A_1, A_2]
    return 1/((1/500) + sum([1/(j+sum([(X[i]-A[i][j-1])**6 for i in range(0,len(X))])) for j in range(1,26)]))

f5 = TestFn("Shekel's Foxholes", shekel, dimension=2)

# Easom function in 2 dimensions

EASOM = TestFn("Easom function", lambda X: -math.cos(X[0])*math.cos(X[1])*math.exp(-1*(X[0]-math.pi)**2 - (X[1]-math.pi)**2), dimension=2)