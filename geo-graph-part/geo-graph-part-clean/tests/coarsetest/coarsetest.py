# -*- coding: utf-8 -*-
"""
Created on Thu May 31 13:58:33 2018

@author: A671118
"""
import sys
from pathlib import Path
sys.path.append(str(Path('...').absolute().parent))
import igraph as ig
import coarse
import random

g = ig.Graph()
nData = 8
nTasks = 4
nVertices = nData + nTasks
maxClusters = nData+nTasks

for i in range(nData):
    g.add_vertex("d%s" % i)

for i in range(nTasks):
    g.add_vertex("t%s" % i)
    
g.vs["volume"] = [random.random() if i < nData else 0. 
    for i in range(nVertices)]
g.vs["index"] = [i for i in range(nVertices)]

g.add_edges([('d0', 't0'), ('d0', 't1'), ('d1', 't0'),
             ('d1', 't2'), ('d2', 't2'), ('d3', 't0'),
             ('d3', 't2'), ('d4', 't1'), ('d4', 't3'),
             ('d6', 't3'), ('d7', 't3'), ('d4', 't0'),
             ('t0', 'd6'), ('t1', 'd5'), ('t2', 'd6'),
             ('t3', 'd5')])

nEdges = len(g.es)
g.es["weight"] = [g.vs["volume"][e.source] + g.vs["volume"][e.target] 
    for e in g.es]

passed = True
for nClusters in range(1, maxClusters + 1):
    g2 = coarse.coarsed_graphs(g, nClusters, quad_assign=True)
    print("Test %s" % nClusters)
    print(len(g2[-1].vs), nClusters, "%s steps" % len(g2))
    print("Passed: %s\n" % (len(g2[-1].vs) == nClusters))
    if len(g2[-1].vs) != nClusters:
        passed=False
print("All tests passed: %s" % passed)