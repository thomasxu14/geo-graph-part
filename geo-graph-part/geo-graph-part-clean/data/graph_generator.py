# -*- coding: utf-8 -*-
"""
Created on Wed May 30 14:29:30 2018

@author: A671118
"""

import igraph as ig
import random
import numpy as np

def generate_graph(nData, nTasks, edgeProb):
    g = ig.Graph()
    nVertices = nTasks + nData
    for i in range(nData):
        g.add_vertex("d%s" % i)
    
    for i in range(nTasks):
        g.add_vertex("t%s" % i)
        
    g.vs["volume"] = [round(random.random(), 2) if i < nData else 0. 
        for i in range(nVertices)]
    
    for t in range(nTasks):
        outputData = random.randint(0, nData-1)
        g.add_edge("d{0}".format(outputData), "t{0}".format(t), weight 
                           = g.vs[nData + t]["volume"] + g.vs[outputData]["volume"])
        for d in range(nData):
            if d != outputData and random.random() <= edgeProb:
                g.add_edge("d{0}".format(d), "t{0}".format(t), weight 
                           = g.vs[nData + t]["volume"] + g.vs[d]["volume"])
    return g

def generate_barabasi(nVertices, outgoingEdges):
    g = ig.Graph.Barabasi(nVertices, outgoingEdges)
        
    g.vs["volume"] = [round(random.random(), 2) for i in range(nVertices)]
    g.vs["name"] = ["v%s" % v.index for v in g.vs]
    g.es["weight"] = [g.vs["volume"][e.source] + g.vs["volume"][e.target] 
        for e in g.es]
    return g

def generate_realistic(nData, nTasks, avgDataDep, stddDataDep):
    g = ig.Graph()
    nVertices = nTasks + nData
    for i in range(nData):
        g.add_vertex("d%s" % i)
    
    for i in range(nTasks):
        g.add_vertex("t%s" % i)
    g.vs["volume"] = [max(round(random.random(), 2), 0.1) if i < nData else 0. 
        for i in range(nVertices)]
    
    for t in range(nTasks):
        nDataDep = min(max(2, int(round(random.gauss(avgDataDep, stddDataDep)))), nData)
        dataDep = random.sample(list(range(nData)), nDataDep)
        for d in dataDep:
            g.add_edge("d{0}".format(d), "t{0}".format(t), weight 
                           = g.vs[nData + t]["volume"] + g.vs[d]["volume"])
    return g