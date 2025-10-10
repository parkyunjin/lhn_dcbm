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
    G1 = nx.read_graphml("/result/dcsbm.graphml")
    G2 = nx.read_graphml("/result/balanced_sbm.graphml")

    # Calculate curvature
    G1 = balancedformancurv(G1)
    G2 = balancedformancurv(G2)

    # Save updated graphs back to GraphML
    nx.write_graphml(G1, "/result/dcsbm_curved.graphml")
    nx.write_graphml(G2, "/result/balanced_sbm_curved.graphml")