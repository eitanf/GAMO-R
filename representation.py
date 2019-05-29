"""
Author: Hrishee Shastri
May 2019


Implementation of an interface between the lower level representations and real numbers
"""
from sympy.combinatorics.graycode import GrayCode 
import math
import random

class Representation:
    """
    Takes a representation function r that maps from the a set of b-bit bitstrings
    to some real interval.
    """
    def __init__(self, repFn, name):
        self._rep = repFn   # bitstr maps to number
        self._invRep = {v: k for k, v in repFn.items()} # number maps to bitstr
        self._name = name 

    def to_num(self, bitstr):
        return self._rep[bitstr]

    def to_bitstr(self, num):
        return self._invRep[num]

    def num_bits(self):
        return len(next(iter(self._rep)))

    def get_random_bitstr(self):
        return random.choice(list(self._rep))

    def get_name(self):
        return self._name

    def is_valid(self, i):
        # Checks if a bitstring i is valid in the real interval. If i is a number,
        # checks if i has a valid bit representation
        return (i in self._rep) or (i in self._invRep)

    def __str__(self):
        return str(self._rep)



def initializeEncodings(encoding, interval):
    """
    Creates the representation function r between an encoding scheme and the real interval.

    encoding -- an ordered list of base 2 bitstrings 
    interval -- a 3-tuple specifying (start, end, step) inclusive, e.g. (-5, 5, 0.1)

    returns a dictionary representing r
    """
    if not isValidInterval(interval):
        raise ValueError("bad interval")

    start = interval[0]
    end = interval[1]
    step = interval[2]

    #number of decimal places is the max number of dp used in either start, end, or step
    dp = max(str(start)[::-1].find('.'),str(end)[::-1].find('.'),str(step)[::-1].find('.'), 0)  

    if len(encoding) < abs((end - start)/step)+1:
        raise ValueError("More items in the interval than there are bitstrings in encoding")

    if start > end: 
        end, start = start, end
        step *= -1

    rep = {}
    j = round(start, dp)
    i = 0

    while j <= end:
        rep[encoding[i]] = j
        i+=1
        j = round(j+step, dp)

    return rep 



def isValidInterval(interval):
    """
    checks if a real interval is valid
    """
    if len(interval) != 3:
        return False 
    start = interval[0]
    end = interval[1]
    step = interval[2]
    return (start < end and step  > 0) or (start > end and step < 0) 



def numBitsToEncodeInterval(interval):
    """
    returns the minimum number of bits b needed to encode all items in a given real interval
    """
    if not isValidInterval(interval):
        raise ValueError("bad interval")
    start = interval[0]
    end = interval[1]
    step = interval[2]

    size = abs((end - start)/step)+1
    return math.ceil(math.log(size, 2))



def generateGrayRepresentation(interval):
    """
    returns gray code as an instance of the Representation class 
    for a given real interval to be used in optimization
    """
    b = numBitsToEncodeInterval(interval)
    gc = list(GrayCode(b).generate_gray())
    grayRep = initializeEncodings(gc, interval)
    return Representation(grayRep, "gray")



def generateBinaryRepresentation(interval):
    """
    returns binary code as an instance of the Representation class 
    for a given real interval to be used in optimization
    """
    b = numBitsToEncodeInterval(interval)
    bc = []
    for i in range(0,2**b):
        binstr = bin(i)[2:]
        bc.append(('0'*(b-len(binstr))+binstr))
    binRep = initializeEncodings(bc, interval)
    return Representation(binRep, "binary")



def generateCustomRepresentation(interval):
    """
    If you want, this is where you can define and implement your own  
    encoding schemes and test their GA performance on the test functions 
    """
    pass  

# br = generateBinaryRepresentation((0,1023,1))
# print(br.num_bits())