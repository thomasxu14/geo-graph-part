# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 14:28:22 2018

@author: A671118
"""

import random
import numpy as np
import igraph as ig
# import erdosrenyi as gdat
import data.testgraph3 as gdat
from data.graph_generator import *
from LP import partitionLP

random.seed(0)

nClusters = 11
clusters = range(nClusters)
tolerance = 0.2

costexp = 30
coststdd = 10

cost = [[0 for c1 in range(nClusters)] for c2 in range(nClusters)]
for c1 in range(nClusters):
    for c2 in range(c1+1, nClusters):
        cost[c1][c2] = max(round(random.gauss(costexp, coststdd),2),10)
        cost[c2][c1] = cost[c1][c2]
cost = np.array(cost)

# cost = np.array([[0., 1.], [1.,0.]])

print("Cost matrix:\n%s" % cost)

nData = 8
nTasks = 3
edgeProb = 0.2
outgoingEdges = 2
avgDataDep = 2
stddDataDep = 0.5

print("Generating graph...")
#g = gdat.g
# g = generate_barabasi(nData + nTasks, outgoingEdges)
g = generate_realistic(nData, nTasks, avgDataDep, stddDataDep)
# g = ig.Graph.Read_Pickle("data/graph_example")
# ig.plot(g)
print("End of generation.")

color_dict_vertex = {0: "blue", 1: "red", 2: "green", 3: "pink", 4: "orange" }
placement = partitionLP(cost, g, tolerance, relaxed=0, quad_assign=True)
partition = []
for vertex in range(len(g.vs)):
    maxProb = 0
    maxPlacement = -1
    for c in range(nClusters):
        if placement[vertex][c] > maxProb:
            maxProb = placement[vertex][c]
            maxPlacement = c
    partition.append(maxPlacement)