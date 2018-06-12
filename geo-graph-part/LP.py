# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 15:57:35 2018

@author: A671118
"""

import igraph as ig
from ortools.linear_solver import pywraplp

def partitionLP(cost, nClusters, g, tolerance, color_dict_vertex = None, relaxed=0):
    clusters = range(nClusters)
    nVertices = len(g.vs)
    vertices = range(nVertices)
    nEdges = len(g.es)
    edges = [(e.source, e.target) for e in g.es]
    volumes = g.vs["volume"]
    totalVolume = sum(volumes)
    loadPerCluster = totalVolume * (1 + tolerance) / nClusters
    edgeWeight = {(i,j): g.es[g.get_eid(i,j)]["weight"] for (i,j) in edges}
    
    # Instantiate a mixed-integer solver, naming it SolveIntegerProblem.
    solver = pywraplp.Solver('SolveIntegerProblem',
                           pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
    if relaxed == 1:
        edgecut = {(i,j,k,l): solver.NumVar(0, +solver.infinity(), 'e_{0}_{1}_{2}_{3}'.format(i,j,k,l)) for (i,j) in edges for k in clusters for l in clusters}
    else:
        edgecut = {(i,j,k,l): solver.BoolVar('e_{0}_{1}_{2}_{3}'.format(i,j,k,l)) for (i,j) in edges for k in clusters for l in clusters}
    v = [[solver.BoolVar('v_{0}_{1}'.format(i,k)) for k in clusters] for i in vertices]
    
    objective = solver.Objective()
    for (i,j) in edges:
        for k in clusters:
            for l in clusters:
                objective.SetCoefficient(edgecut[(i,j,k,l)], edgeWeight[(i,j)] * cost[k,l])
    objective.SetMinimization()
    
    # All the vertices are in one partition
    distributed = [solver.Constraint(1., 1.) for i in vertices]
    for i in vertices:
        for k in clusters:
            distributed[i].SetCoefficient(v[i][k], 1.)
            
    # Balance out the partition load
    balance = [solver.Constraint(-solver.infinity(), loadPerCluster) for k in clusters]
    for k in clusters:
        for i in vertices:
            balance[k].SetCoefficient(v[i][k], volumes[i])
    
    # Linearisation of v_i_k * v_j_l = e_i_j_k_l
    lin1 = {(i,j): [[solver.Constraint(-solver.infinity(), 0) for l in clusters] for k in clusters] for (i,j) in edges}
    lin2 = {(i,j): [[solver.Constraint(-solver.infinity(), 0) for l in clusters] for k in clusters] for (i,j) in edges}
    lin3 = {(i,j): [[solver.Constraint(-1, solver.infinity()) for l in clusters] for k in clusters] for (i,j) in edges}
    for (i,j) in edges:
        for k in clusters:
            for l in clusters:
                lin1[(i,j)][k][l].SetCoefficient(edgecut[(i,j,k,l)], 1)
                lin1[(i,j)][k][l].SetCoefficient(v[i][k], -1)
                lin2[(i,j)][k][l].SetCoefficient(edgecut[(i,j,k,l)], 1)
                lin2[(i,j)][k][l].SetCoefficient(v[j][l], -1)
                lin3[(i,j)][k][l].SetCoefficient(edgecut[(i,j,k,l)], 1)
                lin3[(i,j)][k][l].SetCoefficient(v[i][k], -1)
                lin3[(i,j)][k][l].SetCoefficient(v[j][l], -1)
                
    result_status = solver.Solve()
    # The problem has an optimal solution.
    assert result_status == pywraplp.Solver.OPTIMAL
    assert(solver.VerifySolution(1e-7, True))
    
    print('Number of variables =', solver.NumVariables())
    print('Number of constraints =', solver.NumConstraints())
    
    # The objective value of the solution.
    print('Optimal objective value = %d' % solver.Objective().Value())
    print()
    print("Placement:")
    placement = []
    
    for i in vertices:
        for k in clusters:
            if v[i][k].solution_value() > 0.5:
                #print("{0} --> {1}".format(g.vs["name"][i],k))
                placement.append(k)
                
    cutEdges = []
    totalCost = 0
    for e in g.es:
        if placement[e.source] != placement[e.target]:
            cutEdges.append(True)
#            print("totalCost: %s" % totalCost)
#            print("edge weight: %s" % g.es[e.index]["weight"])
#            print("cluster costs: %s" % cost[placement[e.source], placement[e.target]])
            totalCost += g.es[e.index]["weight"] * cost[placement[e.source], placement[e.target]]
        else:
            cutEdges.append(False)
    print("Total cost: %s" % totalCost)
                
    visual_style = {}
    if color_dict_vertex is not None:
        visual_style["vertex_color"] = [color_dict_vertex[placement[i]] for i in vertices]
        visual_style["edge_color"] = ["black" if cutEdges[i] else color_dict_vertex[placement[g.es[i].source]] for i in range(nEdges)]
    ig.plot(g, **visual_style)
    return placement, visual_style