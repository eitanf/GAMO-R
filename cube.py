import networkx as nx
import representation as rp
import pickle
import numpy as np
from sympy.combinatorics.graycode import GrayCode 


"""
Generates Non greedy gray code using hamiltonian cycles on the hypercube
"""

def hamilton(G, start):
    """
    from https://gist.github.com/mikkelam/ab7966e7ab1c441f947b
    """
    F = [(G,[start])]
    n = G.number_of_nodes()
    while F:
        graph,path = F.pop()
        confs = []
        for node in graph.neighbors(path[-1]):
            conf_p = path[:]
            conf_p.append(node)
            conf_g = nx.Graph(graph)
            conf_g.remove_node(path[-1])
            confs.append((conf_g,conf_p))
        for g,p in confs:
            if len(p)==n:
                return p
            else:
                F.append((g,p))
    return None



def computeSingleBitLocality(rep):
    """
    compute and returns the single bit locality l_r for a given representation function rep.
    """
    rep = rep.get_rep()
    strSum = 0
    for bitString in list(rep.keys()):
        bitSum = 0
        b = len(bitString)
        for pos in range(b):
            mutatedBS = list(bitString)
            if bitString[pos] == "1":
                mutatedBS[pos] = "0"
            else:
                mutatedBS[pos] = "1"
            mutatedBS = ''.join(mutatedBS)
            bitSum += abs(rep[mutatedBS] - rep[bitString])
        strSum+=bitSum
    return strSum/((2**b)*b)

def generateNonGreedyGray(bits=5):
    assert bits >= 3, "ngg only exists for l >= 3"
    g  = nx.hypercube_graph(bits)
    l = bits - 3
    zero = tuple([0]*l)+(0,0,0)
    one = tuple([0]*l)+(0,0,1)
    two = tuple([0]*l)+(0,1,1)
    three = tuple([0]*l)+(1,1,1)
    g.remove_node(zero)
    g.remove_node(one)
    g.remove_node(two)
    g.remove_node(three)
    cycle = [zero, one, two, three] + hamilton(g, tuple([0]*l) + (1,0,1))

    stringifiedCycle = []

    for node in cycle:
        string = ''
        for coord in node:
            string += str(coord)
        stringifiedCycle.append(string)

    repfn = {stringifiedCycle[i] : i for i in range(len(stringifiedCycle))}
    nongreedy = rp.Representation(repfn, "Non Greedy Gray " + str(bits) + "-bits")
    return nongreedy



def recursive_gray(prev_gray):
    """
    prev_gray : list of b-1 bit gray codes
    """
    L = ["0" + bitstr for bitstr in prev_gray]
    R = prev_gray.copy()
    R.reverse()
    R = ["1" + bitstr for bitstr in R]
    return L + R


def nongreedygray_recursive(bits):
    """
    Generates NGG for high bit numbers, when hamiltonian path search is too slow
    """
    if bits == 3:
        ngg3 = generateNonGreedyGray(3)
        return list(ngg3.get_rep().keys())
    else:
        return recursive_gray(nongreedygray_recursive(bits-1))


# with open("NGG_8.txt", "wb") as f:
#     pickle.dump(nongreedygray_recursive(8), f)

# with open("NGG_10.txt", "wb") as f:
#     pickle.dump(nongreedygray_recursive(10), f)

# with open("NGG_12.txt", "wb") as f:
#     pickle.dump(nongreedygray_recursive(12), f)

# with open("NGG_17.txt", "wb") as f:
#     pickle.dump(nongreedygray_recursive(17), f)





# b = 8
# ngg_rep = generateNonGreedyGray(b)
# l_r = computeSingleBitLocality(ngg_rep)
# print(l_r, ((2**b)-1)/b)
# if l_r > ((2**b)-1)/b:
#     with open("NGG_" + str(b) + ".txt", 'wb') as f:
#         pickle.dump(list(ngg.get_rep().keys()), f)
# [0, 1, 19, 2, 31, 28, 20, 3, 23,26, 24, 25, 22, 27, 21, 4, 13, 14, 18, 15, 30, 29, 17, 16,12, 9, 11, 10,7, 8, 6, 5