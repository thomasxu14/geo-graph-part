# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:10:52 2018

@author: A671118
"""

import igraph as ig
import random
import numpy as np

def coarsen(g, quad_assign=False, vertexLimit=0, group_isolated=False, more_balanced=False):
    ### HEM Method for coarsening
    # Calculate a maximal matching, and returns a new coarseGraph based on the matching
    nVertices = len(g.vs)
    edgeWeight = g.es["weight"]
    # vertexWeight = g.vs["volume"]
    matched = [False for i in range(nVertices)]
    matchedNumber = 0
    matching = []
    vertexCount = nVertices
    if more_balanced:
        vertex_list = np.argsort(g.vs["volume"])
    else:
        vertex_list = random.sample(list(range(nVertices)), nVertices)
    for vertex in vertex_list:
        if quad_assign and vertexLimit > 0 and vertexCount == vertexLimit:
            break
        if not matched[vertex]:
            if group_isolated and not g.vs[vertex].neighbors() and matchedNumber/nVertices <= 0.5:
                randomVertex = random.choice(list(range(nVertices)))
                while matched[randomVertex] or randomVertex == vertex:
                   randomVertex = random.choice(list(range(nVertices)))
                matching.append((vertex, randomVertex))
                matched[vertex] = True
                matched[randomVertex] = True
                matchedNumber += 2
                vertexCount -= 1
            else:
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
                    matchedNumber += 2
                    vertexCount -= 1
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


def coarsed_graphs(g, vertexLimit, quad_assign=False, group_isolated=False):
    graphs = [g]
    if len(g.vs) < vertexLimit:
        return graphs
    while len(graphs[-1].vs) > vertexLimit:
        graphs.append(coarsen(graphs[-1], quad_assign=quad_assign, vertexLimit=vertexLimit, group_isolated=group_isolated, more_balanced=True))
    return graphs