"""
Author: Hrishee Shastri
May 2019

Genetic Algorithm for optimization of scalar functions with vector input. 
"""

from chromosome import *


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
    key -- min for function minimization and max for function maximization 
    """
    
    assert popsize%2 == 0, "popsize is not even"
    assert popsize > 0, "popsize is not positive"
    assert 0 <= mutrate and mutrate <= 1, "invalid mutation rate"
    assert 0 <= crossrate and crossrate <= 1, "invalid crossover rate"
    assert gens > 0, "num of generations not positive"

    print("Initializing...")

    # Initialize representation 
    REP = rep(interval)

    print(key.__name__.upper() + "IMIZING " + str(fn).upper() + " (" + REP.get_name() + ")")

    # Initialize file
    f = open("output" + "\\" + file + ".txt", 'w')
    f.write(key.__name__.upper() + "IMIZING " + str(fn).upper() + " (" + str(gens) + " gens, " + REP.get_name() + ")" + "\n")
    f.write(str(fn).upper() + "\n")
    f.write(str(gens) + "\n")
    f.write(REP.get_name() + "\n")
    f.write("Gen #" + "\t" + key.__name__ + "ima" + "\n")


    # Initialize random population
    curr_gen = 1
    POP = []
    dim = fn.get_input_dimension()

    for i in range(0, popsize):
        vec = ""
        for n in range(dim):
            vec += REP.get_random_bitstr()
        POP.append(Chromosome(REP, vec))

    assert len(POP) == popsize, "POP has incorrect number of elements"


    # evaluate population 
    print("Evolving...")
    FITNESS_MAP = {chrom:chrom.eval_fitness(fn) for chrom in POP}
    start_best = key_with_fittest_val_dict(FITNESS_MAP, key)
    curr_best_sol = start_best.to_real_vec()
    curr_best_fitness = FITNESS_MAP[start_best]

    print("Gen count: " + str(curr_gen) + "\t" + "X = " + str(curr_best_sol) + ", " + "fn(X) = " + str(curr_best_fitness))
    f.write(str(curr_gen) + "\t" + str(curr_best_fitness) + "\n")

    # Evolve
    while curr_gen < gens:
        curr_gen += 1
        child_POP = []
        # binary tournament selection
        t_size = 2           
        for i in range(popsize//2):
            parent1 = tournament_selection(POP, t_size, FITNESS_MAP, key)
            parent2 = tournament_selection(POP, t_size, FITNESS_MAP, key)
            while parent2 == parent1:
                parent2 = tournament_selection(POP, t_size, FITNESS_MAP, key)

            if random.uniform(0,1) < crossrate:
                (child1, child2) = parent1.crossover(parent2)
            else:
                child1, child2 = parent1, parent2

            if random.uniform(0,1) < mutrate:
                child1.mutate()
            if random.uniform(0,1) < mutrate:
                child2.mutate()

            child_POP.append(child1)
            child_POP.append(child2)

        POP = child_POP
        assert len(POP) == popsize, "popsize not maintained after next generation"
        FITNESS_MAP = {chrom:chrom.eval_fitness(fn) for chrom in POP}
        candidate = key_with_fittest_val_dict(FITNESS_MAP, key)
        
        # check for improvement
        if key(curr_best_fitness,FITNESS_MAP[candidate]) == FITNESS_MAP[candidate] and curr_best_fitness != FITNESS_MAP[candidate]:
            curr_best_sol = candidate.to_real_vec()
            curr_best_fitness = FITNESS_MAP[candidate]
            print("Gen count: " + str(curr_gen) + "\t" + "X = " + str(curr_best_sol) + ", " + "fn(X) = " + str(curr_best_fitness))
        f.write(str(curr_gen) + "\t" + str(curr_best_fitness) + "\n")

    # return most fit solution 
    print("All " + str(curr_gen) + " generations completed")

    # input vector that produced most optimal scalar output
    X = curr_best_sol 
    # the optima found by the GA     
    y = curr_best_fitness  

    print(key.__name__ + "imum found: " + str([X,y]))
    f.write(key.__name__ + "imum found: " + str(X) + ", " + str(y) + "\n")
    f.write("\n")
    f.close()
    print()
    return (X,y)




