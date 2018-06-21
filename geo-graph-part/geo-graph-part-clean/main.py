# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 11:28:34 2018

@author: A671118
"""

import igraph as ig
import coarse
import random
from graph_utils import visualize, delete_isolated
import data.graph_generator as gen
import numpy as np
from LP import partitionLP
import refining as ref

nData = 100
nTasks = 30
edgeProb = 0.10
outgoingEdges = 2
avgDataDep = 4
stddDataDep = 2
tolerance = 0.01

g = gen.generate_realistic(nData, nTasks, avgDataDep, stddDataDep)
# g = gen.generate_graph(nData, nTasks, edgeProb)
#g = gen.generate_barabasi(nData + nTasks, outgoingEdges)
# from data.testgraph3 import *

# delete_isolated(g)

volumes = g.vs["volume"]
maxVolume = max(volumes)

nVertices = len(g.vs)
nEdges = len(g.es)

visual_style = {}

original_layout = g.layout("fr")
visual_style["layout"] = original_layout
visual_style["vertex_size"] = [8 * (1 + volumes[i]/maxVolume) for i in range(nVertices)]
visual_style["bbox"] = (500, 500)
visual_style["margin"] = 20
visualize(g, visual_style, original_graph=True)
ig.plot(g, "viz/graph.png", **visual_style)
ig.plot(g, **visual_style)

nClusters = 10
clusters = range(nClusters)
costexp = 30
coststdd = 30

cost = [[0 for c1 in range(nClusters)] for c2 in range(nClusters)]
for c1 in range(nClusters):
    for c2 in range(c1+1, nClusters):
        cost[c1][c2] = max(round(random.gauss(costexp, coststdd),2),2.)
        cost[c2][c1] = cost[c1][c2]
cost = np.array(cost)
print("Cost matrix:\n%s" % cost)

vertexLimit = 30
graphs = coarse.coarsed_graphs(g, nClusters, quad_assign=True, group_isolated=True)

cg = graphs[-1]
print("Coarsest graphs volumes: %s" % cg.vs["volume"])
volumes = cg.vs["volume"]
maxVolume = max(volumes)
loadPerCluster = sum(volumes) * (1+tolerance) / nClusters
loadLimits = [loadPerCluster for i in range(nClusters)]

nVertices = len(cg.vs)
vertices = range(nVertices)
nEdges = len(cg.es)

visual_style = {}
coarse_layout = cg.layout("fr")
visual_style["layout"] = coarse_layout
visual_style["vertex_size"] = [8 * (1 + volumes[i]/maxVolume) for i in range(nVertices)]
visual_style["bbox"] = (500, 500)
visual_style["margin"] = 20
print("The graph has been coarsened %s times." % (len(graphs)-1))
ig.plot(cg, "viz/coarsest_graph.png", **visual_style)
ig.plot(cg, **visual_style)

# cluster_color = [[random.randint(0,255) for i in range(3)] for cl in range(nClusters)]
cluster_color = ["#%06x" % random.randint(0, 0xFFFFFF) for cl in range(nClusters)]
print("Colors: %s" % cluster_color)

placement = partitionLP(cost, cg, tolerance, relaxed=0, quad_assign=True)
partition = []
for vertex in range(len(cg.vs)):
    maxProb = 0
    maxPlacement = -1
    for c in range(nClusters):
        if placement[vertex][c] > maxProb:
            maxProb = placement[vertex][c]
            maxPlacement = c
    partition.append(maxPlacement)

checkCost = ref.edgeCut(cg, partition, cost, nClusters)
print("Verifying total cost: " + str(checkCost))

visual_style["vertex_size"] = 20
visual_style["bbox"] = (600, 600)
visual_style["margin"] = 60

volumes = cg.vs["volume"]
maxVolume = max(volumes)
visual_style["vertex_size"] = [10 + volumes[i]*50/(2*maxVolume) for i in range(len(cg.vs))]
visualize(cg, visual_style, partition=partition, color_dict=cluster_color, colored=True)

ig.plot(cg, "viz/initial_partition.png", **visual_style)
ig.plot(cg, **visual_style)

partitions = [partition]
partCosts = [checkCost]

for i in range(len(graphs)-1, 0, -1):
    fg = graphs[i-1]
    cg = graphs[i]
    newPart = coarse.uncoarsen(cg, fg, partitions[0])
    volumes = fg.vs["volume"]
    maxVolume = max(volumes)
    
    newPart = ref.K_L(fg, newPart, nClusters, cost, loadLimits, more_balanced=True)
    partitions.insert(0, newPart)

    visual_style = {}
    layout = fg.layout("fr")
    visual_style["layout"]=layout
    visual_style["vertex_size"] = [8 * (1 + volumes[i]/maxVolume) for i in range(len(fg.vs))]
    visualize(fg, visual_style, partition=newPart, color_dict=cluster_color, original_graph=(i==1), colored=True)

    ig.plot(fg, "viz/partition%s.png" % i, **visual_style)
    print("Free Space: %s" % ref.freeSpace(fg, newPart, loadLimits))
    partCosts.insert(0, ref.edgeCut(fg, newPart, cost, nClusters))
print("\nFinal partition: %s" % partitions[0])
visual_style["layout"] = original_layout
ig.plot(g, **visual_style)

print("Final cost: %s" % ref.edgeCut(g, partitions[0], cost, nClusters))
print("Partial costs: %s" % partCosts)