import os
import pandas as pd
import cv2
import numpy as np
import networkx as nx
from math import comb

# Load adjacency matrix
adj = pd.read_csv('adjacencyMatrix.csv', header=None)
file_name = r'your_filename_for_the_spatial_graph'
G = nx.from_pandas_adjacency(adj)
n = len(adj)

# Initialize local and global metric dictionaries
localMetrics = {}
globalMetrics = {}

# Basic network parameters
cellIndex = list(G.nodes())
localMetrics['Object Index'] = np.array(cellIndex) # object identifier

num_edges = G.number_of_edges()
globalMetrics['Total Edges'] = num_edges

num_nodes = G.number_of_nodes()
globalMetrics['Total Nodes'] = num_nodes

kn = np.sum(adj, axis=1)  # vector with degrees of each node
aveK = np.mean(kn)  # average degree
globalMetrics['Average Degree'] = aveK

# Degree-related metrics

# Network density of graph = normalized average degree
networkDensity = 2*num_edges/(num_nodes*(num_nodes-1))
globalMetrics['Network Density'] = networkDensity

# Variance of degree sequence normalized by maximum possible degree (n - 1)
norm_kn = kn/(num_nodes - 1)
varK = np.var(norm_kn)
globalMetrics['Variance in Degree'] = varK

# Average of normalized neighbor degree sequence
nds = nx.average_neighbor_degree(G)
neighDegreeSequence = list(nds.values())
localMetrics['Average Neighbor Degree'] = np.array(neighDegreeSequence)

normNeighDegreeSequence = np.array(neighDegreeSequence)/(num_nodes-1)
avgNeighborDegree = np.mean(normNeighDegreeSequence)
globalMetrics['Average Neighbor Degree'] = avgNeighborDegree

# Variance of normalized neighbor degree sequence
varNeighborDegree = np.var(normNeighDegreeSequence)
globalMetrics['Variance in Neighbor Degree'] = varNeighborDegree

# Network heterogeneity reflects tendency of hub nodes
networkHeterogeneity = np.std(kn)/aveK
globalMetrics['Network Heterogeneity'] = networkHeterogeneity

# Degrees
def degrees(adj):
    indeg = adj.sum(axis=1)
    outdeg = (adj.T).sum(axis=1)
    g = nx.from_pandas_adjacency(adj)

    if nx.is_directed(g) == True:
        deg = indeg + outdeg
    else: deg = indeg + np.diag(adj).T

    return [deg,indeg,outdeg]

localMetrics['Degrees'] = np.array(degrees(adj)[0])

### k-neighbors ###
# (neighborhood of nodes that are 'k' links away from a given node)
# http://strategic.mit.edu/docs/matlab_networks/kneighbors.m
def kneighbors(adj, ind, k):
    adjk = np.array(adj, copy=True)
    for i in range(k - 1):
        adjk = np.multiply(adjk, adj)
    kneigh = np.nonzero(adjk[ind, :])[0]  # returns indices of k-neighborhood
    return kneigh

# Subgraphs
# http://strategic.mit.edu/docs/matlab_networks/subgraph.m
# alternatively: nx.subgraph
# https://networkx.org/documentation/stable/reference/generated/networkx.classes.function.subgraph.html#networkx.classes.function.subgraph
def subgraph(adj, S):
    # adj_sub = np.array(adj)[np.ix_(S, S)]
    adj_sub = np.array(adj)[S, S]
    return adj_sub


# Clustering coefficient

# From Matlab:
# function [C1,C2, C] = clustCoeff(A)
# %CLUSTCOEFF Compute two clustering coefficients, based on triangle motifs count and local clustering
# % C1 = number of triangle loops / number of connected triples
# % C2 = the average local clustering, where Ci = (number of triangles connected to i) / (number of triples centered on i)
# % Ref: M. E. J. Newman, "The structure and function of complex networks"
# % Note: Valid for directed and undirected graphs
# https://en.wikipedia.org/wiki/Clustering_coefficient

C = []  # initialize clustering coefficient
ind = np.ones(n)
deg = degrees(adj)[0]  # calculate degrees only once
for i in range(n):
    if deg[i] <= 1:
        C.append(0)
    else: C.append(nx.clustering(G,i))
localMetrics['Clustering coefficient'] = C
globalMetrics['Clustering coefficient (normalized)'] = np.mean(C) # includes deg[i] = 0 unlike Matlab version
# alternatively, globalMetrics['Clustering coefficient (normalized)'] = nx.average_clustering(C)

# Local efficiency (dependent on n, kneighbors, subgraphs)
# a measure of how well information is exchanged within the immediate neighborhood of a node
# The code uses the concept of local efficiency based on the inverse of the shortest path lengths within the neighborhood
# of each node.
# https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.efficiency_measures.local_efficiency.html#networkx.algorithms.efficiency_measures.local_efficiency

LE = [] # initialize local efficiency
LE_truncated = []
def local_efficiency(adj, G):
    for i in range(n):
        neigh = kneighbors(adj,i,1)
        sub_neigh = subgraph(adj,neigh)
        n1 = n2 = len(sub_neigh)
        spl = np.empty([n1, n2])
        for i in range(n1):
            for j in range(n2):
                try:
                    shortest_path = nx.shortest_path_length(G, source=i, target=j)
                    spl[i][j] = (shortest_path)
                except nx.NetworkXNoPath:
                    spl[i][j] = 0
        EfficiencyAllPaths = 1.0 / spl
        EfficiencyAllPaths = EfficiencyAllPaths[np.isfinite(EfficiencyAllPaths)]
        LE.append(np.mean(EfficiencyAllPaths))
        # Find out how to place exceptions for RuntimeWarnings
        # RuntimeWarning: divide by zero encountered in divide
        # RuntimeWarning: Mean of empty slice.
        # RuntimeWarning: invalid value encountered in scalar divide: ret = ret.dtype.type(ret / rcount)
    return LE

G_trunc = nx.from_pandas_adjacency(adj[ind==1]) # only considering nodes with degree > 1
globalMetrics['Local Efficiency'] = nx.local_efficiency(G_trunc) # LE_truncated - localEfficiency variable in CalculateNetworkMetrics.m
localMetrics['Local Efficiency'] = local_efficiency(adj,G)

# Assortativity - Pearson correlation coefficient
# https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.assortativity.degree_pearson_correlation_coefficient.html
assortativity = nx.degree_assortativity_coefficient(G)
globalMetrics['Assortativity'] = assortativity

# Degree distribution
kmax = max(kn) # node with max degree
binCenters = np.arange(0, kmax+1)
degreeCounts, _ = np.histogram(kn, bins=binCenters)

# Rich-club metric for threshold degrees from 1 to maximum degree (n-1)
def rich_club_metric(G,k):
    Nk = np.where(deg>=k)[0] # nodes with degree > k
    if len(Nk) == 0:
        phi = 0
        return phi
    elif len(Nk) == 1:
        return 1 ### NOT SURE cause (Nk-1) = 0, https://www.nature.com/articles/nphys209, https://en.wikipedia.org/wiki/Rich-club_coefficient
    else:
        # adjk = subgraph(adj,Nk)
        adjk = G.subgraph(Nk)
        phi = 2*adjk.number_of_edges()/(len(Nk)*(len(Nk)-1))
        return phi

RCM = np.zeros(kmax-1)
for i in range(1,kmax):
    RCM[i - 1] = rich_club_metric(G, i)

# Average rich-club metric
avgeRichClubMetric = np.mean(RCM)
globalMetrics['Average Rich Club Metric'] = avgeRichClubMetric

# Variance in rich-club metric
varRichClubMetric = np.var(RCM)
globalMetrics['Variance in Rich Club Metric'] = varRichClubMetric

# Motif- and module-related metrics

# Number of isolated nodes
if len(degreeCounts) >= 1:
    nSingleNodes = degreeCounts[1]
else: nSingleNodes = 0
globalMetrics['Isolated Node Count'] = nSingleNodes

# Number of independent pairs of nodes in graph
if len(degreeCounts) >= 2:
    nPairNodes = 0.5*degreeCounts[2]
else: nPairNodes = 0
globalMetrics['Pair Node Count'] = nPairNodes

# Number of triangular loops
if num_nodes <= 2:
    nTriangularLoops = 0
else: nTriangularLoops = np.trace(np.linalg.matrix_power(adj, 3)) / 6 #check linalg function
globalMetrics['Triangular Loop Count'] = nTriangularLoops

# Star motifs (dependent on "deg")
def num_star_motifs(adj,k):
    num = 0
    for i in range(1,len(deg)):
        if deg[i]>=(k-1):
            num = num + comb(deg[i],k-1)
    return num

# Number of 4-node star motifs:
if num_nodes <= 3:
    nStar4 = 0
else: nStar4 = num_star_motifs(adj,4)
globalMetrics['4-star Motif Count'] = nStar4

# Number of 5-node star motif:
if num_nodes <= 4:
    nStar5 = 0
else: nStar5 = num_star_motifs(adj,5)
globalMetrics['5-star Motif Count'] = nStar5

# Number of 6-node star motif:
if num_nodes <= 5:
    nStar6 = 0
else: nStar6 = num_star_motifs(adj,6)
globalMetrics['6-star Motif Count'] = nStar6

# Number of connected components excluding single nodes
#[S, C] = graphconncomp(G) finds the strongly connected components of the graph represented by matrix G using Tarjan's algorithm. A strongly connected component is a maximal group of nodes that are mutually reachable without violating the edge directions. Input G is an N-by-N adjacency matrix that represents a graph. Nonzero entries in matrix G indicate the presence of an edge.
#The number of components found is returned in S (ncc), and C (ccLabels) is a vector indicating to which component each node belongs.

# NetworkX doc: https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.strongly_connected_components.html#networkx.algorithms.components.strongly_connected_components

ncc = [len(c) for c in nx.connected_components(G) if len(c)>1] # returns array with number of nodes in each disconnected subgraph
globalMetrics['Number of Connected Components'] = len(ncc) # nConnectedComponents - number of independent subgraphs

# Average component size
globalMetrics['Average Component Size (normalized)'] = np.mean(ncc) #avgeComponentSize
# Note: doesn't calculate UNIQUE component sizes unlike Matlab version

# Variance in component size
globalMetrics['Variance in Component Size'] = np.var(ncc) # varComponentSize

# Network diameter
# print(nx.diameter(G))
# globalMetrics['Network Diameter'] = nx.diameter(G) # networkDiameter
########## networkx.exception.NetworkXError: Found infinite path length because the graph is not connected #######
# https://stackoverflow.com/questions/33114746/why-does-networkx-say-my-directed-graph-is-disconnected-when-finding-diameter

# Global efficiency
globalMetrics['Network Efficiency (normalized)'] = nx.global_efficiency(G) # globalEfficiency

# Centrality metrics
# not Dangalchev, but Wasserman
# https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.closeness_centrality.html
# https://www.sciencedirect.com/science/article/pii/S0378437105012768
def graphshortestpath(G,i):
    dist_dict = nx.single_source_shortest_path_length(G,i) # first entry is [i,0] (i is the node)
    return dist_dict #np.array(list(dist_dict.items()))

Cl = np.zeros(n)
for i in range(n):
    d = graphshortestpath(G, i)
    d.pop(i) # Remove the distance from node i to itself
    Cl[i] = sum(1 / 2 ** value for value in d.values())/(n-1) # normalized for (n-1)

localMetrics['Closeness Centrality'] = Cl

# Node betweeness
# Number of shortest paths that pass through a node
# https://stackoverflow.com/questions/36552135/python-how-to-compute-the-number-of-shortest-paths-passing-through-one-node
localMetrics['Betweenness Centrality'] = list(nx.betweenness_centrality(G).values())

# Create dataframes and save as .csv
df_LM = pd.DataFrame(localMetrics)
df_LM.to_csv(f'LocalMetrics_{file_name}.csv', index=False) # make unique with file name
df_GM = pd.DataFrame(globalMetrics, index = [0])
df_GM.to_csv(f'GlobalMetrics_{file_name}.csv') # make unique with file name