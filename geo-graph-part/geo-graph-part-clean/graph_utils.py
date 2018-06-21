# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 14:41:26 2018

@author: A671118
"""

def visualize(g, visual_style, partition=None, color_dict=None, original_graph=False, colored=False):
    if colored and color_dict is not None and partition is not None:
        cutEdges = []
        for edge in g.es:
            if partition[edge.source] != partition[edge.target]:
                cutEdges.append(True)
            else:
                cutEdges.append(False)
        visual_style["vertex_color"] = [color_dict[partition[vertex]] for vertex in range(len(g.vs))]
        visual_style["edge_color"] = ["black" if cutEdges[edge] else color_dict[partition[g.es[edge].source]] for edge in range(len(g.es))]
    if original_graph:
        visual_style["vertex_shape"] = ["square" if g.vs[vertex]["name"][0] == 't' else "circle" for vertex in range(len(g.vs))]
    else:
        visual_style["vertex_shape"] = ["circle" for vertex in range(len(g.vs))]
    return

def delete_isolated(g):
    toDelete = []
    for v in g.vs:
        if not v.neighbors():
            toDelete.append(v)
    g.delete_vertices(toDelete)