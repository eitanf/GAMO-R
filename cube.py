import networkx as nx
import representation as rp


"""
Generates Non greedy gray code using hamiltonian cycles on the hypercube
"""

def hamilton(G):
    """
    from https://gist.github.com/mikkelam/ab7966e7ab1c441f947b
    """
    F = [(G,[list(G.nodes())[0]])]
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


def generateNonGreedyGray(bits=5):
    g  = nx.hypercube_graph(bits)
    cycle = hamilton(g)

    stringifiedCycle = []

    for node in cycle:
        string = ''
        for coord in node:
            string += str(coord)
        stringifiedCycle.append(string)

    repfn = {stringifiedCycle[i] : i for i in range(len(stringifiedCycle))}
    nongreedy = rp.Representation(repfn, "Non Greedy Gray " + str(bits) + "-bits")
    return nongreedy





