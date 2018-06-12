# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 10:03:58 2018

@author: A671118
"""

import numpy as np


def K_L(g, partition, nClusters, cost, loadLimits):
    initialCut = edgeCut(g, partition, cost, nClusters)
    bestGainMatrix = allGain(g, partition, nClusters, cost)
    bestPartition = partition[:]
    bestCut = initialCut 
    bestAvailSpace = freeSpace(g, partition, loadLimits)
    
    foundBetterPart = True
    while foundBetterPart:
        foundBetterPart = False
        # Loop goes on as long as the computed currentPartition has improved our solution
        # We set foundBetterPart to False and change it to True if we have found a better partition

        currentPartition = bestPartition[:]
        currentCut = bestCut
        gainMatrix = bestGainMatrix[:]
        availSpace = bestAvailSpace[:]

        movedList = [False for v in range(len(g.vs))]

        foundMove = True
        while foundMove:
            foundMove = False
            # The loop goes on as long as we have found a move during the loop.

            # We compute the sorted list of the indices of the gainMatrix
            sortedGainIndex = [(index // len(g.vs), index % len(g.vs)) for index in np.argsort(gainMatrix, axis=None)]
            nextMove = None
            while nextMove is None and sortedGainIndex:
                # We look for the next suitable move in decreasing order for the gain matrix. This move must not be in
                # movedList already.
                potentialMove = sortedGainIndex.pop()
                while movedList[potentialMove[1]] and sortedGainIndex:
                    potentialMove = sortedGainIndex.pop()
                if movedList[potentialMove[1]]:
                    break

                targetOfMove, vertexToMove = potentialMove[0], potentialMove[1]
                sourceOfMove = currentPartition[vertexToMove]

                # We check if there is a cluster that is overloaded
                overload = False
                for c in range(nClusters):
                    if availSpace[c] < 0:
                        overload = True

                if overload:
                    # If there is an overload, we want the move to not overload further the clusters
                    if availSpace[sourceOfMove] < 0 and availSpace[targetOfMove] - g.vs[vertexToMove]["volume"] > availSpace[sourceOfMove]:
                        nextMove = potentialMove
                else:
                    nextMove = potentialMove
            if nextMove is not None:
                if gainMatrix[nextMove] < 0 and foundBetterPart:
                    break
                foundMove = True
                movedList[vertexToMove] = True
                currentPartition[vertexToMove] = targetOfMove
                availSpace[sourceOfMove] += g.vs[vertexToMove]["volume"]
                availSpace[targetOfMove] -= g.vs[vertexToMove]["volume"]
                currentCut -= gainMatrix[nextMove]
                updateGainMatrix(g, currentPartition, nextMove, gainMatrix, nClusters, cost, sourceOfMove)
                if currentCut < bestCut:
                    if foundBetterPart:
                        bestAvailSpace[sourceOfMove] += g.vs[vertexToMove]["volume"]
                        bestAvailSpace[targetOfMove] -= g.vs[vertexToMove]["volume"]
                        bestPartition[vertexToMove] = targetOfMove
                        updateGainMatrix(g, bestPartition, nextMove, bestGainMatrix, nClusters, cost, sourceOfMove)
                    else:
                        bestPartition = currentPartition[:]
                        bestAvailSpace = availSpace[:]
                        bestGainMatrix = gainMatrix.copy()
                    bestCut = currentCut
                    foundBetterPart = True
    return bestPartition


def gain(g, partition, target, vertex, cost):
    source = partition[vertex]
    totalGain = 0
    for neighbor in g.vs[vertex].neighbors():
        totalGain += g.es[g.get_eid(vertex, neighbor.index)]["weight"] * (
                    cost[source, partition[neighbor.index]] - cost[partition[neighbor.index], target])
    return totalGain


def allGain(g, partition, nClusters, cost):
    gainMatrix = np.zeros((nClusters, len(g.vs)))
    for target in range(nClusters):
        for vertex in range(len(g.vs)):
            gainMatrix[target, vertex] = gain(g, partition, target, vertex, cost)
    return gainMatrix


def freeSpace(g, partition, loadLimits):
    availSpace= loadLimits[:]
    overload = False
    for v in range(len(g.vs)):
        availSpace[partition[v]] -= g.vs[v]["volume"]
        if not overload and availSpace[partition[v]] < 0:
            overload = True
    return availSpace


def edgeCut(g, partition, cost, nClusters):
    cut = 0
    for edge in g.es:
        if partition[edge.source] != partition[edge.target]:
            cut += cost[partition[edge.source], partition[edge.target]] * edge["weight"]
    return cut

def updateGainMatrix(g, partition, nextMove, gainMatrix, nClusters, cost, sourceOfMove):
    vertexToMove = nextMove[1]
    targetOfMove = nextMove[0]
    for neighbor in g.vs[vertexToMove].neighbors():

                    neighIndex = neighbor.index
                    neighPlacement = partition[neighIndex]
#                    print("Gain: %s" % gainMatrix[nextMove])
                    for clTarget in range(nClusters):
#                        print(clTarget, vertexToMove, targetOfMove, neighPlacement, sourceOfMove)
                        edgeWeight = g.es[g.get_eid(vertexToMove, neighIndex)]["weight"]
                        gainMatrix[clTarget, vertexToMove] += edgeWeight * (cost[targetOfMove, neighPlacement] - cost[sourceOfMove, neighPlacement])
                        gainMatrix[clTarget, neighIndex] += edgeWeight * (cost[neighPlacement, targetOfMove] - cost[neighPlacement, sourceOfMove]
                                                                          - cost[clTarget, targetOfMove] + cost[clTarget, sourceOfMove])
    return