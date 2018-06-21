# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 14:28:22 2018

@author: A671118
"""

import sys
from pathlib import Path
sys.path.append(str(Path('...').absolute().parent))
import random
import numpy as np
import igraph as ig
import data.testgraph3 as gdat
import data.graph_generator as gen
from LP import partitionLP

random.seed(0)

nClusters = 3
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

print("Cost matrix:\n%s" % cost)

nData = 22
nTasks = 10
edgeProb = 0.2
outgoingEdges = 2
avgDataDep = 3
stddDataDep = 0.5

print("Generating graph...")
#g = gdat.g
# g = gen.generate_barabasi(nData + nTasks, outgoingEdges)
g = gen.generate_realistic(nData, nTasks, avgDataDep, stddDataDep)
# g = ig.Graph.Read_Pickle("data/graph_example")
ig.plot(g)
print("End of generation.")

color_dict_vertex = {0: "blue", 1: "red", 2: "green", 3: "pink", 4: "orange" }
placement = partitionLP(cost, g, tolerance)
partition = []
for vertex in range(len(g.vs)):
    maxProb = 0
    maxPlacement = -1
    for c in range(nClusters):
        if placement[vertex][c] > maxProb:
            maxProb = placement[vertex][c]
            maxPlacement = c
    partition.append(maxPlacement)
print("Placement: %s\n" % placement)
