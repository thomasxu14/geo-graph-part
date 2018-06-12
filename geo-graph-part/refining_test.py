# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 11:00:22 2018

@author: A671118
"""

import igraph as ig
import numpy as np
import refining as ref
import random
from data.graph_generator import *

nData = 60
nTasks = 26
edgeProb = 0.25

outgoingEdges = 2

avgDataDep = 4
stddDataDep = 3

#g = generate_realistic(nData, nTasks, avgDataDep, stddDataDep)
g = generate_barabasi(nData + nTasks, outgoingEdges)

toDelete = []
for v in g.vs:
    if not v.neighbors():
        toDelete.append(v)
g.delete_vertices(toDelete)

volumes = g.vs["volume"]
maxVolume = max(volumes)

nVertices = len(g.vs)
nEdges = len(g.es)

nClusters = 3
clusters = range(nClusters)


costexp = 30
coststdd = 30

cost = [[0 for c1 in range(nClusters)] for c2 in range(nClusters)]
for c1 in range(nClusters):
    for c2 in range(c1+1, nClusters):
        cost[c1][c2] = max(round(random.gauss(costexp, coststdd),2),2.)
        cost[c2][c1] = cost[c1][c2]
cost = np.array(cost)

# cost = np.array([[0., 1.], [1.,0.]])

print("Cost matrix for clusters:\n%s\n" % cost)

randomPart = [random.randint(0, nClusters-1) for i in range(len(g.vs))]
print("Random partition:\n%s\n" % randomPart)

tolerance = 0.2
loadPerCluster = sum(volumes) * (1+tolerance) / nClusters
loadLimits = [loadPerCluster for i in range(nClusters)]

cutValue = 0
for edge in g.es:
    if randomPart[edge.source] != randomPart[edge.target]:
       cutValue += cost[randomPart[edge.source], randomPart[edge.target]] * edge["weight"]
    
print("Value of the cut for the random partition: %s" % cutValue)
print("Checking value:%s" % ref.edgeCut(g, randomPart, cost, nClusters))


gainMatrix = ref.allGain(g, randomPart, nClusters, cost)
#print(gainMatrix)

randomVertex = random.choice(list(range(len(g.vs))))
randomCluster = random.choice(list(range(nClusters)))
sourceOfMove = randomPart[randomVertex]
print("Random move: ({0}, {1})".format(randomVertex, randomCluster))

gain = ref.gain(g, randomPart, randomCluster, randomVertex, cost)
print("Associated gain: %s" % gain)
randomPart[randomVertex] = randomCluster

ref.updateGainMatrix(g, randomPart, (randomCluster, randomVertex), gainMatrix, nClusters, cost, sourceOfMove)

print("New edge cut: %s" % ref.edgeCut(g, randomPart, cost, nClusters))
print("New calculated edge cut: %s" % (cutValue - gain))

checkMatrix = ref.allGain(g, randomPart, nClusters, cost)
print(gainMatrix - checkMatrix)

refinedPart = ref.K_L(g, randomPart, nClusters, cost, loadLimits)
print(ref.edgeCut(g, refinedPart, cost, nClusters))