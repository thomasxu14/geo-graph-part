# -*- coding: utf-8 -*-
"""
Created on Wed May 30 10:42:57 2018

@author: A671118
"""

import igraph as ig

g = ig.Graph()
g.add_vertices(5)
g.add_edges([(0,1), (0,2), (1,3), (2,3), 
             (3,4), (4,2)])
nEdges = len(g.es)
nVertices = len(g.vs)

g.vs["name"] = ["%s" % i for i in range(nVertices)]
g.vs["volume"] = [1,1,1,3,3]
g.vs["index"] = [i for i in range(nVertices)]

g.es["weight"] = [g.vs["volume"][e.source] + g.vs["volume"][e.target] 
    for e in g.es]

g.es["weight"] = [g.vs["volume"][e.source] + g.vs["volume"][e.target] 
    for e in g.es]