# -*- coding: utf-8 -*-
"""
Created on Wed May 30 14:29:30 2018

@author: A671118
"""

import igraph as ig
import random
import numpy as np

g = ig.Graph()
nData = 50
nTasks = 20
nVertices = nData + nTasks
edgeProb = 1/np.sqrt(nData)

for i in range(nData):
    g.add_vertex("d%s" % i)

for i in range(nTasks):
    g.add_vertex("t%s" % i)
    
g.vs["volume"] = [random.random() if i < nData else 0. 
    for i in range(nVertices)]

for t in range(nTasks):
    outputData = random.randint(0, nData-1)
    g.add_edge("d{0}".format(outputData), "t{0}".format(t))
    for d in range(nData):
        if d != outputData and random.random() <= edgeProb:
            g.add_edge("d{0}".format(d), "t{0}".format(t))

nEdges = len(g.es)
g.es["weight"] = [g.vs["volume"][e.source] + g.vs["volume"][e.target] 
    for e in g.es]