"""
Author: Hrishee Shastri
May 2019

Chromosome class: using a specific 2-bit representation for the actual GA. 
A chromosome is a sequence of bitstrings where each bitstring corresponds to
a specific real number in the interval. 

The number of bitstrings in the sequence is equal to the dimension of the function's
input vector. For example, if we were to optimize a function 

f: R^3 --> R

then there would be 3 bitstrings in the chromosome. In other words, a chromosome 
represents a real valued vector.

This file contains much of the genetic operators needed for the GA. 
"""
import random
import numpy

class Chromosome:
    def __init__(self, rep, vector):
        """
        vector -- sequence of bitstrings or a sequence of real numbers (will implicitly convert to bitstrings)
                  note: the sequence of bitstrings will be a string (e.g. "0110101010") where vector entries can be
                  distinguished by knowing the number of bits to encode a real number in the interval. A 
                  sequence of real numbers will be a list of real numbers in the interval, which will be converted 
                  to a bitstring.
        rep -- representation function to be used (Representation object)
        """
        self._b = rep.num_bits()   # number of bits to encode a single real number
        self._rep = rep
        self._vec = ""

        if type(vector) == list or type(vector) == tuple:
            for n in vector:
                self._vec += self._rep.to_bitstr(n)
        else:
            self._vec = vector

    def to_real_vec(self):
        """
        converts bitstring sequence to a real-valued vector
        """
        return [self._rep.to_num(self._vec[i:i+self._b]) for i in range(0, len(self._vec), self._b)]

    def eval_fitness(self, fn):
        """
        computes the fitness of self based on the function fn being optimized. fn is a TestFn object. 
        Note that self.to_real_vec() is the genotype to phenotype mapping, and fn is the phenotype to R mapping.
        """
        return fn.eval(self.to_real_vec())

    def is_valid(self, string=None):
        """
        checks whether self is a valid chromosome by checking if each bitstring is valid.
        string -- if left None, checks if self is a valid chromosome. Otherwise checks if 
                  the string argument could be a valid bitstring sequence. Uses the same 
                  individual bitstring length and representation function as self.
        """
        if string is None:
            string = self._vec

        for i in range(0, len(string), self._b):
            if not self._rep.is_valid(string[i:i+self._b]):
                return False
        return True

    def crossover(self, partner):
        """
        Returns two child chromosomes created from self and a partner chromosome.
        The technique used here is one point crossover.
        """
        p1, p2 = str(self), str(partner)
        assert(len(p1) == len(p2))
        point = random.randint(0,len(p1))
        child1 = p1[:point] + p2[point:]
        child2 = p2[:point] + p1[point:]
        assert(len(child1) == len(p1))
        assert(len(child2) == len(p1))
        return [Chromosome(self._rep, child1), Chromosome(self._rep, child2)]


    def mutate(self, pm):
        """
        multi-bit mutation. Called after mutation rate check is made.
        returns new mutated chromosome. pm = mutation rate
        """
        flip = lambda b : "0" if b == "1" else "1"
        l = len(self._vec)
        mutindiv = []
        for i in range(l):
            if random.uniform(0,1) <= pm:
                mutindiv.append(flip(self._vec[i]))
            else:
                mutindiv.append(self._vec[i])
#        if ''.join(mutindiv) == self._vec:
#            return self
        return Chromosome(self._rep, ''.join(mutindiv))



    def performance_value(self, fmap, f_prime, key):
        """
        u(x). f_max is determined by scaling window.
        """
        if key == min:
            return f_prime - fmap[self]
        else:
            return fmap[self] - f_prime

    def copy(self):
        """
        returns a copy of itself as a new object
        """
        return Chromosome(self._rep, self._vec)

    def __str__(self):
        return self._vec




# Now we define some functions which will help in the GA
def wheel_selection(pop, fmap, f_prime, key):
    """
    Selects two individuals accordining to a fitness proportion distribution 
    pop -- list of chromosomes 
    fmap - fitness map
    key -- min if minimizing fitness and max if maximizing fitness
    """
    w = [indiv.performance_value(fmap, f_prime, key) for indiv in pop]
    # print(f_prime, [i.eval_fitness(fn) for i in pop])
    s = sum(w)
    if s == 0:
        return numpy.random.choice(pop, 2)
    return numpy.random.choice(pop, 2, p = [i/s for i in w])




def tournament_selection(pop, k, fmap, key):
    """
    Selects the most fit individual out of a random k-subset of pop 
    pop -- list of chromosomes (s
    k -- tournament size. Must be lte than len(pop), greater than 0
    fmap -- fitness map, where chromosomes map to their fitness
    key -- min if minimizing fitness and max if maximizing fitness
    """
    assert (0 < k and k <= len(pop)), "invalid tournament size"

    ksubset = []
    popcopy = pop[:]
    for i in range(k):
        ind = random.randint(0,len(popcopy)-1) 
        chrom = popcopy[ind]
        ksubset.append(chrom)
        popcopy.pop(ind) 

    assert len(ksubset) == k, "tournament size not met"

    most_fit = ksubset[0]
    for chrom in ksubset:
        david = fmap[chrom]
        goliath = fmap[most_fit] 
        if key(david, goliath) == david:
            most_fit = chrom 

    return most_fit 

def key_with_fittest_val_dict(d, key):
    """
    just a general purpose function to grab the dictionary key with most fit value

    d -- the dictionary
    key -- min or max
    """
    opt = None
    for k, v in d.items():
        if opt is None or key(d[k], d[opt]) == d[k]:
            opt = k
    return opt


