"""
Author: Hrishee Shastri
May 2019
Runs the GA on De Jong's (1975) test functions with a number of trials and plots the result. 
The main objective of this project is to compare different encoding schemes (e.g. gray vs binary) in optimization tasks.
Accordingly, this program optimizes a function using GAs and plots the Fitness vs Generation curve for gray and binary on the same plot.
"""

import testFunctions as tf
import representation as rp
import multiprocessing as mp
from pathos.multiprocessing import ProcessingPool as Pool
from optimizationGA import GA_SEARCH

# Global constants
GRAY_CODE = rp.generateGrayRepresentation
BINARY_CODE = rp.generateBinaryRepresentation
CUSTOM_CODE = rp.generateCustomRepresentation  # custom encoding scheme. You can define it in representation.py 
NGG_CODE = rp.generateNGG
UBL_CODE = rp.generateUBL


# Parameters suggested by Grefenstette (1986) for optimization of De Jong's (1975) five function test suite
# https://www.academia.edu/6763441/Optimization_of_Control_Parameters_for_Genetic_Algorithms
# However, there are many other variables to take into account, such as type of crossover, type of representation,
# tournament size, etc. 
m = 0.01 # mutation rate
c = 0.95 # crossover rate
p = 30   # population size

g = 100 # no. of generations. Doesnt actually do anything because we run the GA until 5000 fitness evals.

NUM_RUNS = 1000
# minimization
key = min
def main():
    # 5 trials
    # writes online performance to text file function_representation_trial#.txt

    # endpoint in interval must be end - step to make sure the number of discrete points on each axis is a power of 2
            # e.g. is literature says the search space is -5.12 <= x <= 5.12 with resolution \delta x = 0.01, input
            #       (-5.12, 5.11, 0.01) as the interval


    pool = Pool(mp.cpu_count())
    jobs = []

    funcs = [tf.f1, tf.f2, tf.f3, tf.f4, tf.f5]
    ranges = [(-5.12,5.11,0.01), (-2.048,2.047,0.001), (-5.12,5.11,0.01), (-1.28, 1.27, 0.01), (-65.536, 65.535, 0.001)]
    search = lambda i,j : GA_SEARCH(m, c, p, g, GRAY_CODE, "f" + str(j) + "_BRG_T" + str(i), funcs[j-1], ranges[j-1], min)

    for j in range(1, len(funcs)+1):
        print(str(funcs[j-1]))
        for i in range(1,NUM_RUNS+1):
            job = pool.apipe(GA_SEARCH, m, c, p, g, NGG_CODE, "f" + str(j) + "_NGG_T" + str(i), funcs[j-1], ranges[j-1], min)
            jobs.append(job)
            job = pool.apipe(GA_SEARCH, m, c, p, g, UBL_CODE, "f" + str(j) + "_UBL_T" + str(i), funcs[j-1], ranges[j-1], min)
            jobs.append(job)

            job = pool.apipe(GA_SEARCH, m, c, p, g, GRAY_CODE, "f" + str(j) + "_BRG_T" + str(i), funcs[j-1], ranges[j-1], min)
            jobs.append(job)
            job = pool.apipe(GA_SEARCH, m, c, p, g, BINARY_CODE, "f" + str(j) + "_BIN_T" + str(i), funcs[j-1], ranges[j-1], min)
            jobs.append(job)

    for job in jobs:
        job.get()
#    pool.close()
#    pool.join()

if __name__ == "__main__":
    main()
