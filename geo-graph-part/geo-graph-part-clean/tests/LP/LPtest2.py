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
from LP import partitionLP, partitionLP_homogeneous

random.seed(0)

nClusters = 10
clusters = range(nClusters)
tolerance = 0.2

nData = 7
nTasks = 3
edgeProb = 0.2
outgoingEdges = 2
avgDataDep = 2
stddDataDep = 0.5

print("Generating graph...")
#g = gdat.g
# g = generate_barabasi(nData + nTasks, outgoingEdges)
g = generate_realistic(nData, nTasks, avgDataDep, stddDataDep)
# ig.plot(g)
print("End of generation.")

color_dict_vertex = {0: "blue", 1: "red", 2: "green", 3: "pink", 4: "orange" }
placement = partitionLP_homogeneous(g, tolerance, nClusters)
partition = []
for vertex in range(len(g.vs)):
    maxProb = 0
    maxPlacement = -1
    for c in range(nClusters):
        if placement[vertex][c] > maxProb:
            maxProb = placement[vertex][c]
            maxPlacement = c
    partition.append(maxPlacement)