"""
Various functions to compute and compare locality metrics for bitstring to integer representations

Todo: turn these into methods for the Representation class 
"""
from sympy.combinatorics.graycode import GrayCode 
from scipy.special import comb
import math
import random
from representation import *
import cube 

def phenClosed(b):
    n = 2**b
    return (1/6)*(n-1)*(n)*(n+1)

def genClosed(b):
    return b*(2**(2*(b-1)))

def hamming(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def computeDistanceDistortion(rep):
    """
    computes distance distortion of a representation, as defined in Rothlauf 2nd ed. 
    page 84. 
    - Our phenotypic distance metric is simply euclidean distance.
    - Our genotypic distance metric is hamming distance
    """
    dd = 0
    rep = rep.get_rep()
    keys = list(rep.keys())

    for i in range(len(keys)):
        innersum = 0   
        for j in range(i+1, len(keys)):
            d_p = abs(rep[keys[i]] - rep[keys[j]])
            d_g = hamming(keys[i], keys[j])
            #print(str(d_p) + " - " + str(d_g) + " = " + str(abs(d_p - d_g)))
            innersum += abs(d_p - d_g)

        dd += innersum

    return dd/comb(len(rep), 2, exact=True)

def computeDistanceDistortionTriangle(rep):
    """
    computes distance distortion of a representation, as defined in Rothlauf 2nd ed. 
    page 84. 
    - Our phenotypic distance metric is simply euclidean distance.
    - Our genotypic distance metric is hamming distance
    """
    dd = 0
    rep = rep.get_rep()
    keys = list(rep.keys())
    phen = 0
    gen = 0

    for i in range(len(keys)):
        innersum = 0   
        for j in range(i+1, len(keys)):
            d_p = abs(rep[keys[i]] - rep[keys[j]])
            d_g = hamming(keys[i], keys[j])
            #print(str(d_p) + " - " + str(d_g) + " = " + str(abs(d_p - d_g)))
            innersum += d_p - d_g
            phen += d_p
            gen += d_g 

        dd += innersum

    return dd/comb(len(rep), 2, exact=True)



def s_cardinality(rep):
    sprime = 0
    sdoubleprime = 0
    negsum = 0

    rep = rep.get_rep()
    keys = list(rep.keys())
    for i in range(len(keys)):
        for j in range(i+1, len(keys)):
            d_p = abs(rep[keys[i]] - rep[keys[j]])
            d_g = hamming(keys[i], keys[j])
            if d_p >= d_g:
                sprime += 1
            else:
                sdoubleprime += 1
                negsum += d_g - d_p

    return (negsum, sdoubleprime)




def computeDistanceDistortionNoGenotype(rep):
    """
    computes distance distortion of a representation, as defined in Rothlauf 2nd ed. 
    page 84. 
    - Our phenotypic distance metric is simply euclidean distance.
    - Our genotypic distance metric is hamming distance
    """
    dd = 0
    rep = rep.get_rep()
    keys = list(rep.keys())

    for i in range(len(keys)):
        innersum = 0   
        for j in range(i+1, len(keys)):
            d_p = abs(rep[keys[i]] - rep[keys[j]])
            #d_g = hamming(keys[i], keys[j])
            #print(str(d_p) + " - " + str(d_g) + " = " + str(abs(d_p - d_g)))
            innersum += abs(d_p)

        dd += innersum

    return dd/comb(len(rep), 2, exact=True)



def computeRothlaufLocality_1sted(rep):
    """
    Computes Rothlauf's definition for locality of a representation as in the 1st edition
    """
    rep = rep.get_rep()
    keys = list(rep.keys())
    dp_min = 1
    dg_min = 1

    dm = 0
    test = 0
    for i in range(len(keys)):
        for j in range(len(keys)):
            if abs(rep[keys[i]] - rep[keys[j]]) == dp_min:
                dm += abs(hamming(keys[i], keys[j]) - dg_min)
    return dm/2 #divide by 2 since we counted each pair twice


def computeRothlaufLocality_2nded(rep):
    """
    Computes Rothlauf's definition for locality of a representation as in the 1st edition
    """
    rep = rep.get_rep()
    keys = list(rep.keys())
    dp_min = 1
    dg_min = 1

    dm = 0
    for i in range(len(keys)):

        for j in range(len(keys)):
            if abs(hamming(keys[i], keys[j])) == dg_min:
                dm += abs(abs(rep[keys[i]] - rep[keys[j]]) - dp_min) 

    return dm/2 #divide by 2 since we counted each pair twice



def computeSingleBitLocality(rep):
    """
    compute and returns the single bit locality l_r for a given representation function rep.
    """
    rep = rep.get_rep()
    strSum = 0
    for bitString in list(rep.keys()):
        bitSum = 0
        b = len(bitString)
        for pos in range(b):
            mutatedBS = list(bitString)
            if bitString[pos] == "1":
                mutatedBS[pos] = "0"
            else:
                mutatedBS[pos] = "1"
            mutatedBS = ''.join(mutatedBS)
            bitSum += abs(rep[mutatedBS] - rep[bitString])
        strSum+=bitSum
    return strSum/((2**b)*b)










