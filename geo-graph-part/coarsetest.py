# -*- coding: utf-8 -*-
"""
Created on Thu May 31 13:58:33 2018

@author: A671118
"""

import igraph as ig
import coarse
import random

g = ig.Graph()
nData = 8
nTasks = 4
nVertices = nData + nTasks

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

g2 = coarse.coarsen(g)