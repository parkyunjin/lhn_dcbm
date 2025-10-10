## CALCULATE BFC for real data ##
import networkx as nx
from scipy.io import loadmat
import sys
sys.path.insert(1, "/code/curvatures/prev_curv")
from balanced_forman_curvature import BalancedForman

def balancedcurv(G: nx.Graph):
    bfc = BalancedForman(G)
    G = bfc.compute_balancedformancurv()
    return G


if __name__ == '__main__':
    mat = loadmat('/data/caltech.mat')
    adj = mat['A']
    G = nx.from_numpy_array(adj)
    G = balancedcurv(G)
    print("BFC DONE")
    nx.write_gml(G, "/data/caltech_curvs.gml")
    print("Graph saved")

    mat = loadmat('/data/simmons.mat')
    adj = mat['A']
    G = nx.from_numpy_array(adj)
    G = balancedcurv(G)
    print("BFC DONE")
    nx.write_gml(G, "/data/simmons_curvs.gml")
    print("Graph saved")

    mat = loadmat('/data/polblog.mat')
    adj = mat['A']
    G = nx.from_numpy_array(adj)
    G = balancedcurv(G)
    print("BFC DONE")
    nx.write_gml(G, "/data/polblog_curvs.gml")
    print("Graph saved")