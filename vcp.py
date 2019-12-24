"""
Reproduces Eremeev's mutation based tournament selection vertex cover solving GA

The 1999 results
https://link.springer.com/content/pdf/10.1007%2F10721187.pdf (1999)

The 2018 results:
http://delivery.acm.org/10.1145/3250000/3243516/evco_a_00210.pdf?ip=134.10.30.87&id=3243516&acc=ACTIVE%20SERVICE&key=B63ACEF81C6334F5%2E2BBF49B3B72B7709%2E4D4702B0C3E38B35%2E4D4702B0C3E38B35&__acm__=1576988018_104da1ca1a7a5b221ada93f1068dc044
This is virtually identical to the 1999 experiment vcp.py, except we only check 
1000 independent runs of the algorithm
where for each t only one individual g(t)_1 was checked for optimality. Thus for each t we
have a series of 1000 Bernoulli trials. 
"""

import random
import matplotlib.pyplot as plt
# G_d is a graph of d disconnected triangle subgraphs
# Represent G_d as list of undirected edges in [u,v] notation
# This ensures the edges are ordered

def generate_G(d):
    """
    Generates a graph consisting of d disconnected triangle subgraphs
    using the above representation. 
    """
    graph = []
    for n in range(0, 3*d, 3):
        graph.append([n,n+1])
        graph.append([n,n+2])
        graph.append([n+2,n+1])
    return graph        


class VCP:
    def __init__(self, d):
        """
        d : number of disconnected triangle subgraphs to cover
        """
        self.G = generate_G(d) 
        self.E = 3*d
        self.V = 3*d

    def fitness(self, genotype):
        """
        Returns fitness value of input genotype. 
        V - x(genotype) 
        where x(genotype) is the number of vertices covered by the genotype

        Bigger number is better
        """
        assert(len(genotype) == self.E)
        vertices = set()
        for i in range(len(genotype)):
            vertices.add(self.G[i][int(genotype[i])])

        return self.V - len(vertices)

    def tournament_selection(self, s):
        """
        s tournament selection.
        Chooses s individuals, uniformly randomly, from current population.
        Returns the individual with the best fitness 
        """
        return max(random.sample(self.pop, k=min(s,self.N)), key = self.fitness)

    def mutate(self, indiv, pm):
        """
        returns a new individual that has been mutated
        pm : mutation rate
        """
        flip = lambda b : '0' if b == '1' else '1'
        newb = ''
        for i in range(len(indiv)):
            if random.uniform(0,1) < pm: 
                newb += flip(indiv[i])
            else:
                newb += indiv[i]

        assert(len(newb) == len(indiv))
        return newb

    def proportion_of_opt_sols(self):
        """
        Computes and returns the proportion of optimal genotypes in the current population

        optimal solution is one where every triangle subgraph is covered by exactly two edges  
        """
        count = 0
        for indiv in self.pop:
            #check each triangle individually. Each triangle needs to be optimal for the entire cover to be optimal
            for i in range(0,len(indiv),3):
                flag = True
                triangle = set()
                triangle.add(self.G[i][int(indiv[i])])
                triangle.add(self.G[i+1][int(indiv[i+1])])
                triangle.add(self.G[i+2][int(indiv[i+2])])

                if len(triangle) != 2:
                    flag = False
                    break

            count += int(flag)

        return count/len(self.pop)

    def one_indiv_optimal(self):
        """
        Returns True if the first individual in population is optimal. else false. 
        """
        indiv = self.pop[0]
 
        #check each triangle individually. Each triangle needs to be optimal for the entire cover to be optimal
        for i in range(0,len(indiv),3):
            flag = True
            triangle = set()
            triangle.add(self.G[i][int(indiv[i])])
            triangle.add(self.G[i+1][int(indiv[i+1])])
            triangle.add(self.G[i+2][int(indiv[i+2])])

            if len(triangle) != 2:
                flag = False
                break

        return flag



    def run(self, tmax, s, pm, N, bernoulli = False):
        """
        tmax : number of generations
        s : number of individuals selected to participate in tournament selection (only 1 gets picked)
        pm : mutation rate
        N : population size (N individuals each generation)
        bernoulli : whether to reproduce Eremeev 1999 results (False) or Eremeev 2018 results (True)

        outputs data: where data[i] is the proportion of optimal solutions at generation i

        t           
        """
        self.tmax = tmax
        self.s = s
        self.N = N
        self.pm = pm

        self.data = []

        # initial population is one where every triangle subgraoh is covered redunantly, as determined in Eremeev (1999)
        # This means every node is in the vertex cover, which is a genotype of '011' for each triangle subgraph
        self.pop = ['011'*(self.E//3)]*self.N 

        for t in range(tmax):
            self.newpop = []
            for k in range(N):
                g = self.tournament_selection(s)
                self.newpop.append(self.mutate(g,pm))
            assert(len(self.newpop) == len(self.pop))
            self.pop = self.newpop.copy()
            if bernoulli:
                self.data.append(self.one_indiv_optimal())
            else:
                self.data.append(self.proportion_of_opt_sols())

        return self.data 

def vcp_plot(datas, labels, title, x_axis, y_axis, fname, scatter = False):
    """
    datas -- a list of lists of data
    labels -- a list of labels for each list of data (for the legend)
    title -- string title for graph
    x_axis -- string label for x axis
    y_axis -- string label for y axis
    fname -- file name
    scatter -- whether points should be connected (False) or not (True)
    """
    assert(len(datas) == len(labels))

    for i in range(len(datas)):
        data = datas[i]
        xs = list(range(len(data)))
        
        ys = data
        label = labels[i]

        if not scatter:
            plt.plot(xs, ys, label = label, linestyle = 'solid') 
        else:
            plt.scatter(xs, ys)

    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.locator_params(axis='x', nbins=4)  # 4 ticks on x-axis
    plt.title(title)
    plt.legend()
    filename = "graphs" + "\\" + fname + ".png"
    print("Output graph to " + filename)
    plt.savefig(filename)


def average_over(A, d, tmax, s, pm, N, bernoulli = False):
    """
    Runs the VCP EA A times and returns a list L where L[i] is the average proportion of optimal solutions at the ith iteration)

    Other parameters are as they are above, documented in VCP.run

    it prints the data at the end. So if you don't want to run the algorithm again, just copy the printed out data and paste it somewhere
    to graph later. 
    """
    average_data = [0]*tmax
    for i in range(A):
        problem = VCP(d)
        data = problem.run(tmax, s, pm, N, bernoulli)
        average_data = [average_data[i] + data[i] for i in range(len(data))]
        print(i)
    average_data = list(map(lambda x : x/A, average_data))
    print(average_data)
    return average_data


# Reproduce 1999 results
def vary_pop_size(TRIALS = 100):
    # 100 Trials
    # Slows down considerably as N gets larger. 
    data1 = average_over(A = TRIALS, d = 8, tmax = 100, s = 2, pm = 0.01, N = 1)
    data2 = average_over(A = TRIALS, d = 8, tmax = 100, s = 2, pm = 0.01, N = 2)
    data10 = average_over(A = TRIALS, d = 8, tmax = 100, s = 2, pm = 0.01, N = 10)
    data100 = average_over(A = TRIALS, d = 8, tmax = 100, s = 2, pm = 0.01, N = 100)
    vcp_plot([data1,data2,data10,data100], ["N = 1", "N = 2", "N = 10", "N = 100"], 
              title ="(1999) Avg prop of opt sols vs. gen. d = 8, pm = 0.01, s = 2",
              x_axis = 't', 
              y_axis = 'avg. prop. of opt. sols in pop', 
              fname = 'vary_pop_size_Eremeev(1999)')


# Reproduce 1999 results
def vary_tourn_size(TRIALS = 100):
    # 100 Trials
    # Slows down considerably as s gets larger. 
    data1 = average_over(A = TRIALS, d = 6, tmax = 20, s = 1, pm = 0.1, N = 100)
    data2 = average_over(A = TRIALS, d = 6, tmax = 20, s = 2, pm = 0.1, N = 100)
    data10 = average_over(A = TRIALS, d = 6, tmax = 20, s = 10, pm = 0.1, N = 100)
    data50 = average_over(A = TRIALS, d = 6, tmax = 20, s = 50, pm = 0.1, N = 100)
    vcp_plot([data1,data2,data10,data50], ["s = 1", "s = 2", "s = 10", "s = 50"], 
              title ="(1999) Avg prop of opt sols vs. gen. d = 6, pm = 0.1, N = 100",
              x_axis = 't', 
              y_axis = 'avg. prop. of opt. sols in pop', 
              fname = 'vary_tourn_size_Eremeev(1999)')

# Reproduce 2018 Results. 
def vary_pop_size_bernoulli(TRIALS = 1000):
    data1 = average_over(A = TRIALS, d = 8, tmax = 20, s = 2, pm = 0.1, N = 1, bernoulli = True)
    data2 = average_over(A = TRIALS, d = 8, tmax = 20, s = 2, pm = 0.1, N = 2, bernoulli = True)
    data10 = average_over(A = TRIALS, d = 8, tmax = 20, s = 2, pm = 0.1, N = 10, bernoulli = True)
    vcp_plot([data1,data2,data10], ["N = 1", "N = 2", "N = 10"], 
              title ="(2018) Avg prop of opt sols vs. gen. d = 8, pm = 0.1, s = 2",
              x_axis = 't', 
              y_axis = 'avg. prop. of opt. sols in pop (1000 bernoulli)', 
              fname = 'vary_pop_size_Eremeev(2018)')

# Reproduce 2018 Results
def vary_tourn_size_bernoulli(TRIALS = 1000):
    data1 = average_over(A = TRIALS, d = 8, tmax = 20, s = 1, pm = 0.1, N = 100, bernoulli = True)
    data2 = average_over(A = TRIALS, d = 8, tmax = 20, s = 2, pm = 0.1, N = 100, bernoulli = True)
    data10 = average_over(A = TRIALS, d = 8, tmax = 20, s = 10, pm = 0.1, N = 100, bernoulli = True) 

    vcp_plot([data1,data2,data10], ["s = 1", "s = 2", "s = 10"], 
              title ="(2018) Avg prop of opt sols vs. gen. d = 8, pm = 0.1, N = 100",
              x_axis = 't', 
              y_axis = 'avg. prop. of opt. sols in pop (1000 bernoulli)', 
              fname = 'vary_tourn_size_Eremeev(2018)')
















