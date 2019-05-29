"""
Author: Hrishee Shastri
May 2019

Runs the GA on De Jong's (1975) test functions with a number of trials and plots the result. 
The main objective of this project is to compare different encoding schemes (e.g. gray vs binary) in optimization tasks.
Accordingly, this program optimizes a function using GAs and plots the Fitness vs Generation curve for gray and binary on the same plot.
"""

import testFunctions as tf
import representation as rp 
from plot import visualize, average_graph
from optimizationGA import GA_SEARCH

# Global constants
GRAY_CODE = rp.generateGrayRepresentation
BINARY_CODE = rp.generateBinaryRepresentation
CUSTOM_CODE = rp.generateCustomRepresentation  # custom encoding scheme. You can define it in representation.py 

# Number of trials for each function (i.e we run the GA NUM_TRIALS times on a function and take the average performance)
NUM_TRIALS = 30  

# Parameters suggested by Grefenstette (1986) for optimization of De Jong's (1975) five function test suite
# https://www.academia.edu/6763441/Optimization_of_Control_Parameters_for_Genetic_Algorithms
# However, there are many other variables to take into account, such as type of crossover, type of representation,
# tournament size, etc. 
m = 0.05 # mutation rate
c = 0.50 # crossover rate
p = 50   # population size

g = 1000 # no. of generations

# minimization
key = min

def main():
    """
    Runs an experiment with a set number of trials and plots the fitness curve for each function we are optimizing. 
    """

    # Run the experiment
    for i in range(1,NUM_TRIALS+1):
        print("TRIAL NO. " + str(i))

        pg = GA_SEARCH(m, c, p, g, GRAY_CODE, "parabolaG" + str(i), tf.PARABOLA, (-100,100,0.001), key)
        pb = GA_SEARCH(m, c, p, g, BINARY_CODE, "parabolaB"  + str(i), tf.PARABOLA, (-100,100,0.001), key)

        rg = GA_SEARCH(m, c, p, g, GRAY_CODE, "rosenbrockG" + str(i), tf.ROSf, (-5.12, 5.12, 0.01), key)
        rb = GA_SEARCH(m, c, p, g, BINARY_CODE, "rosenbrockB" + str(i), tf.ROSf, (-5.12, 5.12, 0.01), key)

        stg = GA_SEARCH(m, c, p, g, GRAY_CODE, "stepG" + str(i), tf.STEP, (-5.12, 5.12, 0.01), key)
        stb = GA_SEARCH(m, c, p, g, BINARY_CODE, "stepB" + str(i), tf.STEP, (-5.12, 5.12, 0.01), key)

        qg = GA_SEARCH(m, c, p, g, GRAY_CODE, "noisy_quarticG" + str(i), tf.NOISY_QUARTIC, (-1.28, 1.28, 0.01), key)
        qb = GA_SEARCH(m, c, p, g, BINARY_CODE, "noisy_quarticB" + str(i), tf.NOISY_QUARTIC, (-1.28, 1.28, 0.01), key)

        sg = GA_SEARCH(m, c, p, g, GRAY_CODE, "shekelG" + str(i), tf.SHEKEL, (-65.536, 65.536, 0.001),  key)
        sb = GA_SEARCH(m, c, p, g, BINARY_CODE, "shekelB" + str(i), tf.SHEKEL, (-65.536, 65.536, 0.001), key)

        eg = GA_SEARCH(m, c, p, g, GRAY_CODE, "easomG" + str(i), tf.EASOM, (-100, 100, 0.001), key)
        eb = GA_SEARCH(m, c, p, g, BINARY_CODE, "easomB" + str(i), tf.EASOM, (-100, 100, 0.001), key)

    # process data: take averages
    print("Processing data and graphing...")

    pg_g = average_graph(["parabolaG" + str(i) for i in range(1, NUM_TRIALS+1)], "gray")
    pb_g = average_graph(["parabolaB" + str(i) for i in range(1, NUM_TRIALS+1)], "binary") 

    rg_g = average_graph(["rosenbrockG" + str(i) for i in range(1, NUM_TRIALS+1)], "gray") 
    rb_g = average_graph(["rosenbrockB" + str(i) for i in range(1, NUM_TRIALS+1)], "binary") 

    stg_g = average_graph(["stepG" + str(i) for i in range(1, NUM_TRIALS+1)], "gray") 
    stb_g = average_graph(["stepB" + str(i) for i in range(1, NUM_TRIALS+1)], "binary") 

    qg_g = average_graph(["noisy_quarticG" + str(i) for i in range(1, NUM_TRIALS+1)], "gray") 
    qb_g = average_graph(["noisy_quarticB" + str(i) for i in range(1, NUM_TRIALS+1)], "binary")

    sg_g = average_graph(["shekelG" + str(i) for i in range(1, NUM_TRIALS+1)], "gray")  
    sb_g = average_graph(["shekelB" + str(i) for i in range(1, NUM_TRIALS+1)], "binary")

    eg_g = average_graph(["easomG" + str(i) for i in range(1, NUM_TRIALS+1)], "gray")
    eb_g = average_graph(["easomB" + str(i) for i in range(1, NUM_TRIALS+1)], 'binary')


    # Visualization
    visualize([pg_g, pb_g], "3D Parabola minimization (" + str(NUM_TRIALS) + " trials)", 1)
    visualize([rg_g, rb_g], "Rosenbrock Saddle minimization (" + str(NUM_TRIALS) + " trials)", 2)
    visualize([stg_g, stb_g], "Step function minimization (" + str(NUM_TRIALS) + " trials)", 3)
    visualize([qg_g, qb_g], "Noisy Quartic minimization (" + str(NUM_TRIALS) + " trials)", 4)
    visualize([sg_g, sb_g], "Shekel foxholes minimization (" + str(NUM_TRIALS) + " trials)", 5)
    visualize([eg_g, eb_g], "Easom function minimization (" + str(NUM_TRIALS) + " trials)", 6)

    print()
    print("Done")

if __name__ == "__main__":
    main()