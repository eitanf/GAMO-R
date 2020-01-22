"""
A  Markov chain model of the (1+1)-ES experiment.

"""


"""
Long term Markov chain probability
"""

import numpy as np
from representation import *
import math
from scipy.special import comb
import random

def hamming(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def buildTmat(rep, a):
    """
    generates a numpy array (transition matrix) for 
    representtation rep and given a value. The nature of this MC
    is that a is an absorbing state. 



    returns the Transition matrix P and the index g of the global maxima
    Returns [P,g]
    """
    assert(a < 2**rep.num_bits())
    states = list(rep.get_rep().keys())
    xmax = 2**rep.num_bits() -1
    f = [xmax - abs(x - a) for x in range(0,xmax+1)]
    getF = lambda bitstr : f[rep.to_num(bitstr)]

    b = rep.num_bits()

    P = np.zeros((len(states), len(states)))
    for i in range(len(states)):
        state = states[i]
        rowsum = 0
        for j in range(len(states)):
            neighb = states[j]
            h = hamming(state, neighb)
            m = 1/b
            if getF(neighb) > getF(state):
                P[i][j]= (m**h)*((1-m)**(b-h))
                rowsum += P[i][j]
        assert(rowsum <= 1)
        P[i][i] = 1 - rowsum
    return [P,states.index(rep.to_bitstr(a))]      


def meanFPT(P, g, rep, alpha = None, n = 1000):
    """
    Estimates the mean first passage time to state g when starting at a state samples from init distribution alpha
    n is no. of gens
    rep is rep 
    """

    if alpha is None:
        alpha = [1/len(P)]*(len(P))

    FPTs = []
    TRIALS = 10000
    states = list(rep.get_rep().keys())
    for i in range(TRIALS):
        print(i)
        state = random.randint(0,len(states)-1)
        for j in range(n):
            state = np.random.choice(list(range(len(states))), p = list(P[state]))
            if state == g:
                FPTs.append(j+1)
                break        
    m = np.mean(FPTs)
    print(m)
    return m


def theoreticalProbability(P, g, alpha = None, n = 1000):
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

    # long term distribution
    # print("Long term distribution")
    # print(np.dot(np.array(alpha), Pn)) 

    return prob


#####

def main(rep, a):
    P,g = buildTmat(rep,a)
    meanFPT(P, g, rep)




sb = uneitanify(list(range(32)),5, 'binary')

ubl = uneitanify([24, 1, 4, 19, 15, 16, 21, 13, 9, 26, 18, 0, 23, 12, 6, 22, 3, 28, 20, 14, 30, 7, 5, 27, 29, 10, 8, 31, 2, 17, 25, 11], 5, 'worst2')
ngg = uneitanify([0, 1, 19, 2, 31, 28, 20, 3, 23, 26, 24, 25, 22, 27, 21, 4, 13, 14, 18, 15, 30, 29, 17, 16,12, 9, 11, 10, 7, 8, 6, 5],5, 'ngg')

brg = generateGrayRepresentation((0,31,1))


main(ubl, 31)

