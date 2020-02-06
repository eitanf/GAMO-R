# Binary Representations in Genetic and Evolutionary Algorithms

A collection of small programs and genetic and evolutionary algorithms to help understand the effect of locality in binary-integer representations on performance. 

Special python dependencies: sympy, numpy, pathos, matplotlib, pickle, networkx

Special C++ dependencies: Intel TBB Library (only for onemax.cc)

## Overview 
The following lists the most important files for this project. The ones not mentioned here are either simple utilities included by these, or are graphs and output data. 

`representation.py` contains definitions of a Representation object, which is used heavily throughout other Python implementations. It also contains useful functions for initializing common types of representations (e.g. SB, BRG, UBL, NGG), computing various properties (such as no. of local optima), and translating to and from permutation notation. 

`distdistortion.py` contains functions to compute distance distortion and point locality of representations. `locality.cc` is a C++ implementation of the same, for speed.

`cube.py` generates non-greedy Gray codes using Hamiltonian walks on the hypercube. 

`onemax.cc` is the main implementation of the general ONEMAX, for both SA and ES.

`markovAnalysis.py` contains Markov chain models for the Simulated annealing general ONEMAX problem. It computes long term probabilities with and without the temperature parameter.

`ES_markov.py` contains the Markov chain models for the Evolutionary strategies general ONEMAX problem. It estimates first passage times and computes long term probabilities.

`comparison-GA` contains code to reproduce the genetic algorithm from Caruana & Schaffer (1986). All code inside this folder is independent from anything in the top level directory. 

## Reproduction
The paper ('Revisiting Locality in Binary-Integer Representations') contains a high-level overview of the steps to reproduce the experiments in the experimental section. This section details how to run the code. The source files themselves contain more specific documentation, which may answer other questions.
First, download the repository as is, and make sure you have the required dependencies above. For C++ files, we compiled using 

`g++-7 -Wall -Wextra -pedantic -O3 -march=native -std=c++17 [fname].cc -o [fname]`

### Simulated Annealing (SA)
Run `onemax.cc`, changing any parameters in main as desired. Make sure each simulation is evolving as `SA_generation()` (line 329). Data is output each generation to the terminal.
### Evolutionary Strategies (ES)
Run `onemax.cc`, changing any parameters in main as desired. Make sure each simulation is evolving as `ES_generation()` (line 329).
### Genetic Algorithms (GAs)
Run `optimizationGA.py` in the `comparison-GA` folder, changing any parameters as desired. Afterwards, statistics from the runs can be computed using `data_analysis.py` in the same folder. 




