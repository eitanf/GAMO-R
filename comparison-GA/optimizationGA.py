"""
Author: Hrishee Shastri
May 2019
Genetic Algorithm for optimization of scalar functions with vector input. 
"""

from chromosome import *
import os
import math

def GA_SEARCH(mutrate, crossrate, popsize, gens, rep, file, fn, interval, key=min):
    """
    Executes a genetic algorithm to optimize a mathematical function fn. Returns a pair (X,y) where X is an input vector and y is the optimized fn(X)
    mutrate -- mutation rate, between 0 and 1 inclusive
    crossrate -- crossover rate, between 0 and 1 inclusive
    popsize -- positive even integer population size to be maintained throughout iteration
    gens -- a number greater than 0 that specifies the number of generations to iterate through
    rep -- representation function to be used (instance of Representation class). Maps from bitstrings to real numbers in the given interval
           Pass the function object (e.g. GRAY_CODE)
    file -- text file name to write output to (not the same as console output -- file output writes every generation, while
            console output only writes when an improvement has been made)
    fn -- the real valued mathematical function to be optimized, wrapped in a TestFn object. fn : R^n --> R (i.e. vector valued inputs, scalar valued outputs).
    interval -- A 3-tuple (start, end, step) inclusive that constrains the search space for fn. In other words, each entry x_i in the input vector 
                is constrained by x_i \in [start,end] with step increments. Make sure fn is continuous along every point in the interval (e.g. no ZeroDivisionErrors).
    W -- scaling window = 1
    S -- selection strategy = E  
    key -- min for function minimization and max for function maximization 
    """

    assert popsize > 0, "popsize is not positive"
    assert 0 <= mutrate and mutrate <= 1, "invalid mutation rate"
    assert 0 <= crossrate and crossrate <= 1, "invalid crossover rate"
    assert gens > 0, "num of generations not positive"

#    print("Initializing...")

    # Initialize representation 
    REP = rep(interval)

#    print(key.__name__.upper() + "IMIZING " + str(fn).upper() + " (" + REP.get_name() + ")")


    f = open(os.path.join("caruana_data", file + ".txt"), 'w')
    g = open(os.path.join("caruana_data", file + "best_sol" + ".txt"), 'w')

    # Initialize random population
    EVAL_LIMIT = 5000
    EVALS = 0
    curr_gen = 1
    POP = []
    dim = fn.get_input_dimension()

    for i in range(0, popsize):
        vec = ""
        for n in range(dim):
            vec += REP.get_random_bitstr()
        chrom = Chromosome(REP, vec)
        POP.append(chrom)


    assert len(POP) == popsize, "POP has incorrect number of elements"


    # evaluate population 
#    print("Evolving...")
    # Fitness map is not performance value. It is just the evaluation of the objective function to be minimized.
    FITNESS_MAP = {chrom:chrom.eval_fitness(fn) for chrom in POP}

    # scaling window of 1
    if key == min:
        best = math.inf
        f_prime = max(FITNESS_MAP.values())
    else:
        best = -math.inf
        f_prime = min(FITNESS_MAP.values())

    for k in POP:
        # f.write(str(k.performance_value(FITNESS_MAP, f_prime, key)))
        # f.write("\t")
        f.write(str(FITNESS_MAP[k]))
        f.write("\n")
        EVALS += 1

    g.write(str(key(FITNESS_MAP.values())) + "\n")
    # Evolve
    while EVALS < EVAL_LIMIT:
        curr_gen += 1
        child_POP = []
        new_children = []  # new individuals not from previous generation. Child_pop is the entire population that will replace POP.
                            # new_children keeps track of the individuals that are not from previous generation
        for i in range(popsize//2):
            parent1, parent2 = wheel_selection(POP, FITNESS_MAP, f_prime, key)

            if random.uniform(0,1) <= crossrate:
                child1, child2 = parent1.crossover(parent2)
            else:
                child1, child2 = parent1, parent2

            child1 = child1.mutate(mutrate)
            child2 = child2.mutate(mutrate)

            if child1 != parent1 and child1 != parent2:
                new_children.append(child1)
            if child2 != parent1 and child2 != parent2:
                new_children.append(child2)


            child_POP.append(child1)
            child_POP.append(child2)

        # elitist replacement
        best_chrom = key(FITNESS_MAP, key = FITNESS_MAP.get)
        if best_chrom not in child_POP:
            child_POP.append(best_chrom)

        POP = child_POP.copy()

        assert len(POP) == popsize or len(POP) == popsize + 1, "popsize not maintained after next generation"
        FITNESS_MAP = {chrom:chrom.eval_fitness(fn) for chrom in POP}

        # scaling window of 1, so recompute f_prime every generation
        if key == min:
            f_prime = max(FITNESS_MAP.values())
        else:
            f_prime = min(FITNESS_MAP.values())

        for new in new_children:
            # f.write(str(new.performance_value(FITNESS_MAP, f_prime, key)))
            # f.write("\t")
            f.write(str(FITNESS_MAP[new]))
            f.write("\n")
            EVALS += 1
            if EVALS == EVAL_LIMIT:
                break 

        g.write(str(key(FITNESS_MAP.values())) + "\n")

#    print("All " + str(EVALS) + " fitness evals completed")
