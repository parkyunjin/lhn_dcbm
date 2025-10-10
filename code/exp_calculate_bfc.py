import pandas as pd
import networkx as nx
import igraph as ig
import numpy as np

import sys
sys.path.insert(1, "/code/curvatures/prev_curv")
from balanced_forman_curvature import BalancedForman

def balancedformancurv(G: nx.Graph):
    bfc = BalancedForman(G)
    G = bfc.compute_balancedformancurv()
    return G

if __name__ == '__main__':
    
    # Read in graphs
    G1 = nx.read_graphml("/result/graph_n1.graphml")
    G2 = nx.read_graphml("/result/graph_n2.graphml")
    G3 = nx.read_graphml("/result/graph_n3.graphml")
    G4 = nx.read_graphml("/result/graph_n4.graphml")
    G5 = nx.read_graphml("/result/graph_n5.graphml")
    G6 = nx.read_graphml("/result/graph_n6.graphml")
    G7 = nx.read_graphml("/result/graph_n7.graphml")

    # Calculate curvature
    G1 = balancedformancurv(G1)
    G2 = balancedformancurv(G2)
    G3 = balancedformancurv(G3)
    G4 = balancedformancurv(G4)
    G5 = balancedformancurv(G5)
    G6 = balancedformancurv(G6)
    G7 = balancedformancurv(G7)

    # Save updated graphs back to GraphML
    nx.write_graphml(G1, "/result/graph_n1_curved.graphml")
    nx.write_graphml(G2, "/result/graph_n2_curved.graphml")
    nx.write_graphml(G3, "/result/graph_n3_curved.graphml")
    nx.write_graphml(G4, "/result/graph_n4_curved.graphml")
    nx.write_graphml(G5, "/result/graph_n5_curved.graphml")
    nx.write_graphml(G6, "/result/graph_n6_curved.graphml")
    nx.write_graphml(G7, "/result/graph_n7_curved.graphml")