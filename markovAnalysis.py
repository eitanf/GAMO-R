"""
A Theoretical Markov chain model of the Simulated Annealing experiment.
Doesn't simulate the chain -- it computes theoretical long term probabilities of finding the global maxima
with temperature and non-temperature versions of the SA.
"""


"""
Long term Markov chain probability
"""

import numpy as np
from representation import *
import math
import random

def buildTmat(rep, a):
    """
    generates a numpy array (transition matrix) for 
    representtation rep and given a value. The nature of this MC
    is that a is an absorbing state. 

    This is without temperature

    returns the Transition matrix P and the index g of the global maxima
    Returns [P,g]
    """
    assert(a < 2**rep.num_bits())
    states = list(rep.get_rep().keys())
    f = [a - abs(x - a) for x in range(0,2**rep.num_bits())]
    getF = lambda bitstr : f[rep.to_num(bitstr)]

    b = rep.num_bits()

    P = np.zeros((len(states), len(states)))
    for i in range(len(states)):
        state = states[i]
        rowsum = 0
        for neighb in rep.get_neighbors(state):
            j = states.index(neighb)
            if getF(neighb) > getF(state):
                P[i][j]= 1/len(state)
                rowsum += 1/len(state)
        assert(rowsum <= 1)
        P[i][i] = 1 - rowsum
    return [P,states.index(rep.to_bitstr(a))]      



def initTmatTemperature(rep, a, T = 50):
    """
    generates a numpy array (transition matrix) for 
    representtation rep and given a value. The nature of this MC
    is that a is an absorbing state. 

    This is with temperature

    returns the Transition matrix P and the index g of the global maxima
    Returns [P,g]
    """


    assert(a < 2**rep.num_bits())
    states = list(rep.get_rep().keys())
    f = [a - abs(x - a) for x in range(0,2**rep.num_bits())]
    getF = lambda bitstr : f[rep.to_num(bitstr)]

    b = rep.num_bits()

    P = np.zeros((len(states), len(states)))
    for i in range(len(states)):
        state = states[i]
        rowsum = 0
        for neighb in rep.get_neighbors(state):
            j = states.index(neighb)
            if getF(neighb) > getF(state):
                P[i][j]= 1/len(state)
                rowsum += 1/len(state)
            else:
                p = math.exp((getF(neighb)-getF(state))/T)/len(state)
                P[i][j] = p
                rowsum += p
        assert(rowsum <= 1)
        P[i][i] = 1 - rowsum

    return [P,states.index(rep.to_bitstr(a))] 


def theoreticalProbability(P, g, alpha = None, n = 10000000):
    """
    Using MC theory, calculates the proportion of subproblems solved for a given fitness function 
    and representation without temperature. This is for a population size of 1.
    Assumes uniform initial distribution (as in the experiment) and takes large
    matrix powers to approximate long term behaviour. 

    All local maxima are absorbing, so our end goal is to find the long term probability of the chain being 
    absorbed in the global optima and not any of the local optima. 

    Specifically, this function returns

    \alpha \dot P^{n}_{*g} 

    where \alpha is the uniform initial distribution, \dot is dot product, P is the transition matrix, and
    P_{*g} denotes the gth column of P where g is the index coresponding to the global maxima 
    """

    if alpha is None:
        alpha = [1/len(P)]*(len(P))


    Pn = np.linalg.matrix_power(P,n)
    prob = 0
    for i in range(len(P)):
        prob += alpha[i]*Pn[i][g]

    return prob


def theoreticalProbabilityWithTemp(rep, a, alpha = None, n = 10000):
    """
    Here we perform matrix multiplication Pg = P1*P2*...*Pn where
    Pi is the transition matrix for the nonhomogenous markov chain
    at time i. The dot product

    \alpha \dot P^{g}_{*g} 

    is the value that is returned. 
    """
    if alpha is None:
        alpha = [1/2**rep.num_bits()]*(2**rep.num_bits())

    T_0 = 50
    coolRate = 0.995
    currP,g = initTmatTemperature(rep, a, T_0)

    for i in range(n):
        currP = np.matmul(currP, initTmatTemperature(rep,a,T_0*(coolRate)**i)[0])

    prob = 0
    for i in range(len(currP)):
        prob += alpha[i]*currP[i][g]

    return prob
#####

def main(rep, a):
    print("Probabilities...")
    print()
    prob = theoreticalProbabilityWithTemp(rep, a)
    print("With temperature", prob) # with temp
    P,g = buildTmat(rep,a)
    prob = theoreticalProbability(P,g)
    print("Without temperature", prob) #without temp



