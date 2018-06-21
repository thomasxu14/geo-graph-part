# -*- coding: utf-8 -*-
"""
Created on Wed May 30 10:42:57 2018

@author: A671118
"""

import igraph as ig

g = ig.Graph()
g.add_vertices(6)
g.add_edges([(0,1), (0,2), (1,2), (2,3), 
             (3,4), (4,5), (3,5)])
nEdges = len(g.es)
nVertices = len(g.vs)

g.vs["name"] = ["%s" % i for i in range(nVertices)]
g.vs["volume"] = [1 
            for d in range(nVertices)]
g.vs["index"] = [i for i in range(nVertices)]