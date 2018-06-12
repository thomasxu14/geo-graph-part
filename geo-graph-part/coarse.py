# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:10:52 2018

@author: A671118
"""

import igraph as ig
import random

def coarsen(g):
    ### HEM Method for coarsening
    # Calculate a maximal matching, and returns a new coarseGraph based on the matching
    nVertices = len(g.vs)
    edgeWeight = g.es["weight"]
    vertexWeight = g.vs["volume"]
    matched = [False for i in range(nVertices)]
    matching = []
    for vertex in random.sample(list(range(nVertices)), nVertices):
        if not matched[vertex]:
            matchedNeighbor = None
            maxEdgeWeight = 0
            for neighbor in g.vs[vertex].neighbors():
                nindex = neighbor.index
                if not matched[nindex]:
                    currentEdgeWeight = edgeWeight[g.get_eid(vertex, nindex)]
                    if currentEdgeWeight > maxEdgeWeight:
                        matchedNeighbor = nindex
                        maxEdgeWeight = currentEdgeWeight
            if matchedNeighbor != None:
                matching.append((vertex, matchedNeighbor))
                matched[vertex] = True
                matched[matchedNeighbor] = True
    coarseGraph = g.copy()
    coarseGraph.vs["contracted1"] = [v.index for v in g.vs]
    coarseGraph.vs["contracted2"] = [v.index for v in g.vs]
    toDelete = []
    for (v1, v2) in matching:
        toDelete.append(v2)
        coarseGraph.vs[v1]["name"] = coarseGraph.vs[v1]["name"] + coarseGraph.vs[v2]["name"]
        coarseGraph.vs[v1]["contracted1"] = v1
        coarseGraph.vs[v1]["contracted2"] = v2
        coarseGraph.vs[v1]["volume"] = coarseGraph.vs[v1]["volume"] + coarseGraph.vs[v2]["volume"]
        for neighbor in coarseGraph.vs[v2].neighbors():
            if neighbor in coarseGraph.vs[v1].neighbors() and neighbor.index != v1:
                coarseGraph.es[coarseGraph.get_eid(v1, neighbor.index)]["weight"] += coarseGraph.es[coarseGraph.get_eid(v2, neighbor.index)]["weight"]
            elif neighbor.index != v1:
                coarseGraph.add_edge(v1, neighbor.index, weight = 
                                     coarseGraph.es[coarseGraph.get_eid(v2, neighbor.index)]["weight"])
            coarseGraph.delete_edges(coarseGraph.get_eid(v2, neighbor.index))
    coarseGraph.delete_vertices(toDelete)
    return coarseGraph


def uncoarsen(coarseGraph, fineGraph, coarsePartition):
    finePartition = [None for i in range(len(fineGraph.vs))]
    for coarseVertex in coarseGraph.vs:
        cluster = coarsePartition[coarseVertex.index]
        contracted1 = coarseVertex["contracted1"]
        contracted2 = coarseVertex["contracted2"]
        finePartition[contracted1] = cluster
        finePartition[contracted2] = cluster
    return finePartition