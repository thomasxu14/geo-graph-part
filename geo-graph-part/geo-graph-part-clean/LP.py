# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 15:57:35 2018

@author: A671118
"""

from ortools.linear_solver import pywraplp
import time

def partitionLP(cost, g, tolerance, quad_assign=False, relaxed=0):
    nClusters = len(cost)
    clusters = range(nClusters)
    vertices = range(len(g.vs))
    edges = [(e.source, e.target) for e in g.es]
    volumes = g.vs["volume"]
    totalVolume = sum(volumes)
    loadPerCluster = totalVolume * (1 + tolerance) / nClusters
    edgeWeight = {(edge.source, edge.target): edge["weight"] for edge in g.es}
    
    # Instantiate a mixed-integer solver, naming it SolveIntegerProblem.
    print("Initializing solver...")
    solver = pywraplp.Solver('SolveIntegerProblem',
                           pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
    
    # Create the variables v and edgecut
    edgecut = {(i,j,k,l): solver.NumVar(0, 1, 'e_{0}_{1}_{2}_{3}'.format(i,j,k,l)) for (i,j) in edges for k in clusters for l in clusters}
    if relaxed == 1: # Relaxation of int constraint but doesn't yield good results usually.
        v = [[solver.NumVar(0, 1, 'v_{0}_{1}'.format(i,k)) for k in clusters] for i in vertices]
    else:
        v = [[solver.BoolVar('v_{0}_{1}'.format(i,k)) for k in clusters] for i in vertices]
    
    # Set the objective function
    objective = solver.Objective()
    for (i,j) in edges:
        for k in clusters:
            for l in clusters:
                objective.SetCoefficient(edgecut[(i,j,k,l)], edgeWeight[(i,j)] * cost[k,l])
    objective.SetMinimization()
    
    # Constraints
    print("Setting constraints...")
    # All the vertices are in one partition
    distributed = [solver.Constraint(1., 1.) for i in vertices]
    for i in vertices:
        for k in clusters:
            distributed[i].SetCoefficient(v[i][k], 1.)
            
    # Balance out the partition load. If quad_assign is True, then the problem changes
    # and becomes the quadratic assignment problem
    if not quad_assign:
        balance = [solver.Constraint(-solver.infinity(), loadPerCluster) for k in clusters]
        for k in clusters:
            for i in vertices:
                balance[k].SetCoefficient(v[i][k], volumes[i])
    else:
        balance = [solver.Constraint(1, solver.infinity()) for k in clusters]
        for k in clusters:
            for i in vertices:
                balance[k].SetCoefficient(v[i][k], 1)
    
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
                
    print("Solving...")
    start_time = time.time()
    result_status = solver.Solve()
    print("Finished solving in %s seconds." % (time.time() - start_time))
    # The problem has an optimal solution.
    if result_status != pywraplp.Solver.OPTIMAL:
        print("Warning: value not optimal.")
    if not solver.VerifySolution(1e-7, True):
        print("Warning: solution check is False.")
    
    print('Number of integer variables =', sum([len(v[i]) for i in range(len(v))]))
    print('Number of real variables =', len(edgecut))
    print('Number of constraints =', solver.NumConstraints())
    
    # The objective value of the solution.
    print('Optimal objective value = %d' % solver.Objective().Value())
    print()
    # print("Placement:")
    placement = [[v[i][k].solution_value() for k in clusters] for i in vertices]
    return placement

def partitionLP_homogeneous(g, tolerance, nClusters):
    clusters = range(nClusters)
    vertices = range(len(g.vs))
    edges = [(e.source, e.target) for e in g.es]
    volumes = g.vs["volume"]
    totalVolume = sum(volumes)
    loadPerCluster = totalVolume * (1 + tolerance) / nClusters
    edgeWeight = {(edge.source, edge.target): edge["weight"] for edge in g.es}
    
    # Instantiate a mixed-integer solver, naming it SolveIntegerProblem.
    print("Initializing solver...")
    solver = pywraplp.Solver('SolveIntegerProblem',
                           pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
    
    # Create the variables v and edgecut
    contained = {(i,j): solver.NumVar(0, 1, 'contained_{0}_{1}'.format(i,j)) for (i,j) in edges}
    v = [[solver.BoolVar('v_{0}_{1}'.format(i,k)) for k in clusters] for i in vertices]
    
    # Set the objective function
    objective = solver.Objective()
    for (i,j) in edges:
        objective.SetCoefficient(contained[(i,j)], edgeWeight[(i,j)])
    objective.SetMaximization()
    
    # Constraints
    print("Setting constraints...")
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
    
    # Linearisation of contained_i_j = min_k( 1-abs(vik-vjk) )
    lin1 = {(i,j): [solver.Constraint(-solver.infinity(), 1) for k in clusters] for (i,j) in edges}
    lin2 = {(i,j): [solver.Constraint(-solver.infinity(), 1) for k in clusters] for (i,j) in edges}
    for (i,j) in edges:
        for k in clusters:
            lin1[(i,j)][k].SetCoefficient(contained[(i,j)], 1)
            lin1[(i,j)][k].SetCoefficient(v[i][k], -1)
            lin1[(i,j)][k].SetCoefficient(v[j][k], +1)
            lin2[(i,j)][k].SetCoefficient(contained[(i,j)], 1)
            lin2[(i,j)][k].SetCoefficient(v[i][k], +1)
            lin2[(i,j)][k].SetCoefficient(v[j][k], -1)
                
    print("Solving...")
    start_time = time.time()
    result_status = solver.Solve()
    print("Finished solving in %s seconds." % (time.time() - start_time))
    # The problem has an optimal solution.
    if result_status != pywraplp.Solver.OPTIMAL:
        print("Warning: value not optimal.")
    if not solver.VerifySolution(1e-7, True):
        print("Warning: solution check is False.")
    
    print('Number of integer variables =', sum([len(v[i]) for i in range(len(v))]))
    print('Number of real variables =', len(contained))
    print('Number of constraints =', solver.NumConstraints())
    
    # The objective value of the solution.
    print('Optimal objective value = %d' % solver.Objective().Value())
    print()
    # print("Placement:")
    placement = [[v[i][k].solution_value() for k in clusters] for i in vertices]
    return placement