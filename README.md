# Genetic Algorithms for Mathematical Optimization - Representations

A genetic algorithm (GA) for optimizing functions of the form    f: R^n --> R   (i.e. vector input, scalar output). Inspired by Caruana & Schaffer's (1986) study of binary vs. gray code in GAs, this program was conceived with the objective of gaining a better understanding on how different base 2 genotypic representations affect performance. As of now, it supports two different representations (standard binary and binary reflected gray code) but the user can easily define more. 

Executing the main program runs the GA (one with binary genotype, one with gray genotype) on each of De Jong's (1975) test functions a number of times (NUM_TRIALS). The output is printed to the console and to multiple text files (one text file per GA run). The resulting data (average fitness vs. generation #) is averaged and plotted using matplotlib. 

Dependencies: matplotlib, sympy 

#### File overview

`representation.py` -- Abstracts the low level details of the genotype. This is where the user can define and implement their own 
                    representation functions r: {0,1}^b --> R as an instance of a Representation object. 
                    
`chromosome.py` -- Implementation of a chromosome class. Contains much of the genetic operators for the GA. 

`testFunctions.py` -- Implementation of De Jong's test functions for optimization, as well as a couple others. Can easily be augmented if                         desired.

`optimizationGA.py` -- The actual GA. For our genetic operator design choices, we use tournament selection (k=2), single bit mutation, and uniform crossover. 

`plot.py` -- using matplotlib for data visualization. Extracts run data from text files, averages over all trials, and saves graphs to                  /graphs subdirectory.

`main.py` -- main program. Runs the GA with specified parameters on specified functions over a number of trials. Outputs text files to /output and graphs to /graphs

#### Some example graph output
![Fig. 1](/graphs/NoisyQuarticminimization(30trials).png)
![Fig. 2](/graphs/Stepfunctionminimization(30trials).png)
