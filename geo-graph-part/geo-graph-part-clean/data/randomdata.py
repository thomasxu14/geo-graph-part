# -*- coding: utf-8 -*-
"""
Created on Wed May 23 16:35:26 2018

@author: A671118
"""

import random as rand
import numpy.random as nprd
import numpy as np

rand.seed(0)
nprd.seed(1)

nVertices = 10
nClusters = 5
prob = 1/4

volexp = 15
volstdd = 5
costexp = 10
coststdd = 5

volumes = [max(round(rand.gauss(volexp, volstdd),2),0) 
            for d in range(nVertices)]
volumes = np.array(volumes)

cost = [[0 for c1 in range(nClusters)] for c2 in range(nClusters)]
for c1 in range(nClusters):
    for c2 in range(c1+1, nClusters):
        cost[c1][c2] = max(round(rand.gauss(costexp, coststdd),2),0)
        cost[c2][c1] = cost[c1][c2]
cost = np.array(cost)

edges = [(i,j) for i in range(nVertices) for j in range(i+1, nVertices) if rand.random() <= prob]
edges = np.array(edges)
nEdges = len(edges)