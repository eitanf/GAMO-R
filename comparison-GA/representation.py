"""
Author: Hrishee Shastri
May 2019


Implementation of an interface between the lower level representations and real numbers
"""
from sympy.combinatorics.graycode import GrayCode 
import math
import random
import itertools
import pickle


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

    def get_rep(self):
        return self._rep

    def to_bitstr(self, num):
        return self._invRep[num]

    def get_neighbors(self, bitstr):
        """
        returns list of hamming neighbors of bitstr
        """
        flip = lambda b: "1" if b == "0" else "0"
        neighbs = []
        for i in range(len(bitstr)):
            neighbs.append(bitstr[:i] + flip(bitstr[i]) + bitstr[i+1:])
        return neighbs 

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

    assert len(encoding) == round(abs((end - start)/step) + 1), "More items in the interval than there are bitstrings in encoding"

    int(start/step)

    rep = {}
    i = 0
    j = int(start/step)
    for binstr in encoding:
        rep[binstr] = j*step
        j+=1

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

    size = abs((end - start)/step)
    return math.ceil(math.log(size, 2))



def generateGrayRepresentation(interval, b = None):
    """
    returns gray code as an instance of the Representation class 
    for a given real interval to be used in optimization
    """
    if b is None:
        b = numBitsToEncodeInterval(interval)
    gc = list(GrayCode(b).generate_gray())
    grayRep = initializeEncodings(gc, interval)
    return Representation(grayRep, "binary reflected gray")



def generateBinaryRepresentation(interval, b = None):
    """
    returns binary code as an instance of the Representation class 
    for a given real interval to be used in optimization
    """
    if b is None:
        b = numBitsToEncodeInterval(interval)
    bc = []
    for i in range(0,2**b):
        binstr = bin(i)[2:]
        bc.append(('0'*(b-len(binstr))+binstr))

    binRep = initializeEncodings(bc, interval)
    return Representation(binRep, "binary")

def generateUBL(interval, b = None):
    if b is None:
        b = numBitsToEncodeInterval(interval)
    if b == 10:
        fname = "UBL_10.txt"
    elif b == 12:
        fname = "UBL_12.txt"
    elif b == 8:
        fname = "UBL_8.txt"
    elif b == 17:
        fname = "UBL_17.txt"
    else:
        raise ValueError("interval does not support any of the precomputed UBL reps. May need to add")
    with open(fname, 'rb') as f:
        uc = pickle.load(f)

    return Representation(initializeEncodings(uc, interval), "UBL")

def generateNGG(interval, b = None):
    if b is None:
        b = numBitsToEncodeInterval(interval)
    if b == 10:
        fname = "NGG_10.txt"
    elif b == 12:
        fname = "NGG_12.txt"
    elif b == 8:
        fname = "NGG_8.txt"
    elif b == 17:
        fname = "NGG_17.txt"
    else:
        raise ValueError("interval does not support any of the precomputed NGG reps. May need to add")
    with open(fname, 'rb') as f:
        gc = pickle.load(f)

    return Representation(initializeEncodings(gc, interval), "NGG")


def generateModifiedBinaryRepresentation(interval):
    b = numBitsToEncodeInterval(interval)
    bc = []
    for i in range(0,2**b):
        binstr = bin(i)[2:]
        bc.append(('0'*(b-len(binstr))+binstr))
    s1 = random.randrange(0,len(bc) - 1)
    s2 = random.randrange(0,len(bc)-1)
    bc[s1],bc[s2] = bc[s2], bc[s1]
    binRep = initializeEncodings(bc, interval)
    return Representation(binRep, "binary")

def generateRandomRepresentation(interval, name = 'r'):
    """
    returns a random mapping between bitstrings to numbers in the interval
    """
    b = numBitsToEncodeInterval(interval)
    c = list(GrayCode(b).generate_gray())
    random.shuffle(c)
    randRep = initializeEncodings(c, interval)
    return Representation(randRep, name)


def generateWorstRepresentation(nbits, name = 'w'):
    """
    returns an encoding on nbits that has the worst locality using Harpers algorithm.
    """
    b = nbits
    c = list(GrayCode(b).generate_gray())
    rep = {}
    startstr = random.choice(c)
    parity = startstr.count("1") % 2
    rep[startstr] = 0

    sameParity = list(filter(lambda x: x.count("1")%2 == parity, c))
    sameParity.remove(startstr)
    random.shuffle(sameParity)

    oppParity = list(filter(lambda x: x.count("1")%2 != parity, c))
    random.shuffle(oppParity)

    assert(len(sameParity) + 1 + len(oppParity) == len(c))

    for i in range(len(sameParity)):
        rep[sameParity[i]] = i + 1

    for i in range(len(oppParity)):
        rep[oppParity[i]] = i + 2**(b-1)

    rp = Representation(rep, name)

    return Representation(rep, name)
    

def generateAllReps(numbits):
    """
    returns a list of all representations on numbits bits
    """
    c = list(GrayCode(numbits).generate_gray())
    reps = []
    for perm in list(itertools.permutations(c)):
        reps.append(Representation(initializeEncodings(perm, (0, (2**numbits)-1, 1)), 'name'))
    return reps



def generateCustomRepresentation(interval):
    """
    If you want, this is where you can define and implement your own  
    encoding schemes and test their GA performance on the test functions 
    """
    pass  


def countOptima(perm, key = min):
    """
    Counts number of optima (by key = max or min) in a permutation with neighborhood = 2 
    """
    opt = 0
    for p in range(len(perm)):
        if key(perm[p], perm[(p+1) % len(perm)]) == perm[p] and key(perm[p], perm[(p-1) % len(perm)]) == perm[p]:
            opt+=1
    return opt


def allOptimaBitstring(perm, rep, key = max):
    """
    Returns list of local optima -- induced optima (min or max) -- given a function perm and the bitstring representation (neighborhood = len of bitstring)
    perm is the function inducing optima in the bitstrings, as a list
    rep is a representation obj
    """

    flip = lambda x: '1' if x == '0' else '0'

    optima = 0
    rmap = rep.get_rep()
    optlist = []
  
    for b in list(rmap.keys()): 
        # neighbors, including itself
        neighbs = [b[:i] + flip(b[i]) + b[i+1:] for i in range(len(b))]
        opts = [perm[rep.to_num(nb)] <= perm[rep.to_num(b)] for nb in neighbs]
        if all(opts):
            optlist.append(b)
        
    return optlist

def countOptimaBitstring(perm, rep, key=max):
    """
    Counts the number of induced optima (min or max) given a function perm and the bitstring representation (neighborhood = len of bitstring)
    perm is the function inducing optima in the bitstrings, as a list
    rep is a representation obj
    """
    return len(allOptimaBitstring(perm,rep,key))


def optimaFitMetric(a, rep, key = max):
    """
    Metric that computes average fitness difference between local optima and their neighbors.
    a = a value
    rep = rep object
    """
    perm =  [a - abs(x - a) for x in range(0,2**rep.num_bits())]
    optlist = allOptimaBitstring(perm, rep, key)
    globalopt = rep.to_bitstr(a)
    optlist.remove(globalopt)
    if optlist == []:
        return 0
    s = 0
    
    for i in range(len(optlist)):
        for neighb in rep.get_neighbors(optlist[i]):
            s += abs(perm[rep.to_num(optlist[i])] - perm[rep.to_num(neighb)])

    for neighb in rep.get_neighbors(globalopt):
        s -= abs(a - perm[rep.to_num(neighb)])
    return s/(len(optlist)*len(globalopt))




def findOneMaxA(rep1, rep2, b):
    """
    Finds all a values such that the induced number of maxima in rep1 is less than the 
    induced number of maxima in rep2
    """
    avals = []
    perms = [[(2**(b)-1) - abs(x - a) for x in range(2**b)] for a in range(2**b)]
    a = 0
    for perm in perms: 
        b1 = countOptimaBitstring(perm, rep1, max)
        b2 = countOptimaBitstring(perm, rep2, max)
        print(b1, b2)
        if b1 <= b2:
            avals.append(a)
        a += 1

    return avals

def eitanify(rep):
    """
    Rewrites a representation in Eitan's notation
    The  format is a list of integers, where l[i] is the number that the bitstring (when written in binary, represent i) maps to
    """
    binrep = generateBinaryRepresentation((0,2**rep.num_bits()  - 1, 1)).get_rep().keys()
    rep = rep.get_rep()
    return [rep[b] for b in binrep]

def uneitanify(rep, name = ''):
    """
    Rewrites a rep in dict (rep object) notation from Eitan's notation
    """
    d = {}
    lenb = int(math.log(len(rep), 2))
    for i in range(len(rep)):
        x = bin(i)[2:]
        x = '0'*(lenb - len(x)) + x
        d[x] = rep[i]
    return Representation(d, name)

# print(generateWorstRepresentation(17))

# with open("UBL_8.txt", 'rb') as f:
#     l = pickle.load(f)
#     print(l)

