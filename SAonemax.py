"""
Name: Hrishee Shastri
File: SAonemax.py

Reproduces Rothlauf's Simulated Annealing Mutation based search experiment with
binary, reflected gray ('good' gray), and non-greedy gray ('bad' gray).

Experiment appears on page 134, section 5.4.2 of Rothlauf's "Representations For Genetic and Evolutionary Algorithms" 2nd ed.
"""

"""
- Concatenate 10 integer general one-max problems. Each problem has length l. Goal is to get each subproblem to optimal solution 'a'. 
- Fitness of an individual is the sum of fitness of each of the 10 subproblems, calculated using formula 5.2
- Phenotypic individual encoded in gray or binary
- Start temperature T_0 is 50, T_{n+1} = 0.995*T_{n} 
- Mutation flips a random bit
- If offspring better, offspring replaces parent 
- If offspring worse, offspring replaces parent with probability P(T), a function of temperature T -- Boltzmann selection
- For l = 5: a = 15 and a = 31. 
- Perform 100 runs. Each run stopped after 2000 mutation steps 
- Average across the runs and graph the results for each a
"""

import math
import random

import representation as rp
import plot as plt 
import cube as cb
import distdistortion as dd


class SAexperiment:
    def __init__(self, CONCAT, LENGTH, OPTIMUM, STARTTEMP, REPRESENTATION):
        """
        CONCAT = num. of gen-one max problems to concatenate
        LENGTH = length of each subproblem
        OPTIMUM = the solution to the entire problem, a
        STARTTEMP = start temperature
        REPRESENTATION = representation object
        """
        self.concat = CONCAT
        self.L = LENGTH
        self.A = OPTIMUM
        self.T = STARTTEMP
        self.rep = REPRESENTATION

        assert(self.A <= (2**self.L)-1)
        randbit = lambda : random.choice(["0","1"])
        self.x = ''.join([randbit() for i in range(LENGTH*CONCAT)]) # current individual
        self.history = [] # list of individuals, history[0] is first individual, history[n] is nth individual. 
                         # used for post analysis

    def fitness(self, indiv):
        """
        Computes fitness of current individual as the sum of the fitness of each length L subproblem
        Uses eq. 5.2 from the book 

        indiv == individual whose fitness is to be calculated
        """
        xmax = (2**self.L)-1 
        fitnessSum = 0
        count = 0

        for i in range(0, self.L*self.concat, self.L):
            fitnessSum += xmax - abs(self.rep.to_num(indiv[i:i+self.L]) - self.A)
            count += 1

        assert(count == self.concat) # make sure we're summing correctly 

        return fitnessSum 

    def numOfCorrectSubProblems(self, indiv):
        """
        counts number of correct subproblems for a given individual with self.concat problems
        """
        correct = 0
        for i in range(0, self.L*self.concat, self.L):
            correct += int(self.rep.to_num(indiv[i:i+self.L]) == self.A)

        assert(correct <= self.concat)
        return correct

    def mutate(self):
        """
        mutates the current individual by flipping a random bit. 
        Returns the mutated string i.e. the offspring. 
        """
        i = random.choice(range(0,len(self.x)))
        bit = "0" if self.x[i] == "1" else "1"

        offspring = self.x[:i] + bit + self.x[i+1:] 
        assert(len(offspring) == len(self.x))

        return offspring

    def boltzmann(self, offspring):
        """
        Just computes and returns replacement probability when offspring has lower fitness than parent.
        """
        p = math.exp((self.fitness(offspring) - self.fitness(self.x))/self.T)
        #print(p, self.fitness(offspring), self.fitness(self.x))
        assert(p >= 0 and p <= 1)
        return p

    def run(self, steps):
        """
        Runs the SA. 
        steps = number of mutation steps until finished
        Returns the individual at the end of the iteration
        """
        self.history.clear()

        for i in range(steps):
            offspring = self.mutate() 
            if self.fitness(offspring) >= self.fitness(self.x):
                self.x = offspring
            elif random.uniform(0,1) < self.boltzmann(offspring):
                self.x = offspring
            self.T *= 0.995

            self.history.append(self.x)

        return self.x 

    def postProcess(self):
        """
        Produces a list l where l[n] = number of correct sub problems at the nth step.
        Note that history[n] is the individual at the nth step
        """
        l = [self.numOfCorrectSubProblems(self.history[i]) for i in range(len(self.history))]
        graph = plt.Graph(self.rep.get_name())
        for i in range(len(l)):
            graph.add_point((i, l[i]))

        return graph

    def postProcessFitness(self):
        """
        Produces a list l where l[n] = fitness of individuual at the nth step.
        Note that history[n] is the individual at the nth step
        """
        l = [self.fitness(self.history[i]) for i in range(len(self.history))]
        graph = plt.Graph(self.rep.get_name())
        for i in range(len(l)):
            graph.add_point((i, l[i]))

        return graph

def main(OPTIMUM, TRIALS, NUMBITS):
    """
    Runs the SA, and processes the data to produce the aforementioned graphs.
    Graphs are saved to /graphs subdirectory. It will print out the file name 
    and directory when it is done running.

    Runs on BRG, Binary, worst, and NGG for NUMBITS bits.
    """
    L = NUMBITS
    TEMP = 50
    CONCAT = 1
    A = OPTIMUM

    RUNS = TRIALS
    assert(RUNS > 0)

    # Representations 
    binary = rp.generateBinaryRepresentation((0, (2**L)-1, 1))
    reflectedGray = rp.generateGrayRepresentation((0, (2**L)-1, 1))
    nonGreedyGray = cb.generateNonGreedyGray(L)
    worst = rp.uneitanify([9, 26, 20, 13, 22, 14, 8, 21, 17, 3, 5, 19, 11, 25, 31, 6, 16, 4, 1, 30, 15, 23, 24, 7, 2, 27, 29, 0, 18, 12, 10, 28], 'worst, l_r = 16, a = 15, 5 local max')


    bingraphs = []
    reflectedgraygraphs = []
    nongreedygraygraphs = []
    worstgraphs = []

    for i in range(RUNS):
        binexp = SAexperiment(CONCAT, L, A, TEMP, binary)
        reflectedgrayexp = SAexperiment(CONCAT, L, A, TEMP, reflectedGray)
        nongreedygrayexp = SAexperiment(CONCAT, L, A, TEMP, nonGreedyGray)
        worstexp = SAexperiment(CONCAT, L, A, TEMP, worst)

        binexp.run(2000)
        reflectedgrayexp.run(2000)
        nongreedygrayexp.run(2000)
        worstexp.run(2000)

        bingraphs.append(binexp.postProcess())
        reflectedgraygraphs.append(reflectedgrayexp.postProcess())
        nongreedygraygraphs.append(nongreedygrayexp.postProcess())
        worstgraphs.append(worstexp.postProcess())

    # average over RUNS number of runs
    g1 = plt.average_graph(bingraphs, binary.get_name())
    print(g1.get_Ys()[-1])
    g2 = plt.average_graph(reflectedgraygraphs, reflectedGray.get_name())
    print(g2.get_Ys()[-1])
    g3 = plt.average_graph(nongreedygraygraphs, nonGreedyGray.get_name())
    print(g3.get_Ys()[-1])
    g4 = plt.average_graph(worstgraphs, worst.get_name())
    print(g4.get_Ys()[-1])

    plt.plot([g1, g2, g3, g4], 'a = ' + str(A) +' (' + str(RUNS) + ' trials)', 'number of fitness calls', 'number of correct subproblems', 1)

def main_avg_fit(REPS, A_VALS, TRIALS, NUMBITS):
    """
    Runs the SA, and processes the data to produce the aforementioned graphs.
    Graphs are saved to /graphs subdirectory. It will print out the file name 
    and directory when it is done running.

    """
    L = NUMBITS
    TEMP = 50
    CONCAT = 1

    RUNS = TRIALS
    assert(RUNS > 0)

    graphs = []


    for i in range(len(REPS)):
        repg = []
        for j in range(RUNS):
            exp = SAexperiment(CONCAT, L, A_VALS[i], TEMP, REPS[i])

            exp.run(2000)

            repg.append(exp.postProcessFitness())
            print(j)
        graphs.append(repg)

    # save data

    # average over RUNS number of runs
    graphs = [plt.average_graph(graphs[i],REPS[i].get_name()) for i in range(len(REPS))]

    plt.plot(graphs, 'w v b (' + str(RUNS) + ' trials)', 'number of fitness calls', 'average fitness', 1)






# worst = rp.uneitanify([9, 26, 20, 13, 22, 14, 8, 21, 17, 3, 5, 19, 11, 25, 31, 6, 16, 4, 1, 30, 15, 23, 24, 7, 2, 27, 29, 0, 18, 12, 10, 28], 'worst, l_r = 16, a = 15, 5 local max')
# binary = rp.generateBinaryRepresentation((0, (2**5)-1, 1))
# rand = rp.uneitanify([29, 4, 23, 2, 11, 1, 30, 17, 10, 31, 22, 26, 5, 14, 0, 15, 12, 27, 8, 18, 7, 9, 6, 19, 13, 16, 25, 21, 24, 20, 3, 28], 'l_r = 10, a = 15, 5 local max')
# brg = rp.generateGrayRepresentation((0,31,1))
# ngg = cb.generateNonGreedyGray(5)



#main(31,1000,5)
#main_avg_fit([worst, binary], [15,15], 1000, 5)
#main_avg_fit([brg,binary,ngg], [31,31,31], 2000, 5)




