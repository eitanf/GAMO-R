"""
Author: Hrishee Shastri
May 2019

Plotting fitness vs generation number of GAs

Runs with the same optimizing function but different encodings will be shown on the same plot,
to see how different encodings compare
"""

import matplotlib.pyplot as plt 
import statistics as st


class Graph:
    def __init__(self, label=""):
        """
        A class that encapsulates data points for plotting.

        label -- string label for what the graph represents
        """
        self._label = label
        self._x_vals = []
        self._y_vals = []

    def add_point(self, pt):
        self._x_vals.append(pt[0])
        self._y_vals.append(pt[1])

    def get_Xs(self):
        return self._x_vals

    def get_Ys(self):
        return self._y_vals

    def get_label(self):
        return self._label

    def set_Xs(self, xs):
        self._x_vals = xs

    def set_Ys(self, ys):
        self._y_vals = ys

    def __str__(self):
        return str(self._x_vals) + "\n" + str(self._y_vals)

    def __len__(self):
        assert len(self._x_vals) == len(self._y_vals), "X values and Y values do not have same number of elements"
        return len(self._x_vals)



def plot(graphs, title, x_axis, y_axis, fig):
    """
    graphs -- a list of Graph objects to plot on the same plane
    title -- string title for graph
    x_axis -- string label for x axis
    y_axis -- string label for y axis
    fig -- figure number
    """
    plt.figure(fig)
    for g in graphs:
        xs = g.get_Xs()
        ys = g.get_Ys()
        plt.plot(xs, ys, label = g.get_label(), linestyle = 'dashed')

    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.locator_params(axis='x', nbins=4)  # 4 ticks on x-axis
    plt.title(title)
    plt.legend()
    filename = "graphs" + "\\" + ''.join((title.split(' '))) + ".png"
    print("Output graph to " + filename)
    plt.savefig(filename)



def parse(fname):
    """
    Grabs data from text file fname.txt 
    Returns a Graph object containing said data
    """
    f = open("output" + "\\" + fname + '.txt', 'r')
    l = list(f)
    newg = Graph(l[3])
    for line in l[5:]:
        s = line.split('\n')[0]
        s = s.split('\t')
        try:
            newg.add_point((int(s[0]), float(s[1])))
        except ValueError:
            break
    return newg

def average_graph(objs, label):
    """
    objs -- list of filenames or graphs (must all have same length)
    label -- name of resulting graph
    Returns average graph: a graph with average y-value for each x-value across all y-values in each graph 
    Used for experiments where we want to analyze a set number of trials
    """ 
    assert len(objs) > 0, "list contains no graphs"

    if type(objs[0]) != type(Graph("")):
        graphs = [parse(file) for file in objs]
    else:
        graphs = objs

    assert [len(graph) for graph in graphs] == [len(graphs[0])]*len(graphs), "cannot average over graphs with different lengths"

    avg_graph = Graph(label)    
    # x-axis should be same across all graphs (the generation count)
    avg_graph.set_Xs(graphs[0].get_Xs())
    avg_Ys = []
    for i in range(0, len(graphs[0].get_Ys())):
        entry_avg = []
        for graph in graphs:
            entry_avg.append(graph.get_Ys()[i])
        avg_Ys.append(st.mean(entry_avg))
    avg_graph.set_Ys(avg_Ys)
    return avg_graph

def visualize(objs, graphTitle, fig):
    """
    objs -- list of filenames without the .txt suffix OR list of graphs. Graphs the data from each filename or graph on the same plot.
    graphTitle -- string name of plot title
    fig -- figure number
    """
    assert len(objs) > 0, "Cannot visualize nothing"

    if type(objs[0]) != type(Graph("")):
        objs = [parse(file) for file in objs]
    plot(objs, graphTitle, "Generation", "Average Fitness", fig)






