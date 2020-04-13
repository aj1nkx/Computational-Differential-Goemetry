# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 11:29:46 2020

@author: Ajinkya Jayawant
"""

from GraphRicciCurvature.OllivierRicci import OllivierRicci

import networkx as nx
import numpy as np
import math
import importlib
import matplotlib
import matplotlib.pyplot as plt
import logging


def main():
    # grc_example
    jupyter_tutorial()
    return
    
def jupyter_tutorial():
    matplotlib.use('Qt4Agg')
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.ERROR)
    G = nx.karate_club_graph()
    print(nx.info(G))
    orc = OllivierRicci(G, alpha=0.5, verbose="INFO")
    orc.compute_ricci_curvature()
    G_orc = orc.G.copy()  # save an intermediate result
    plt.figure()
    show_results(G_orc)
    # Start a Ricci flow with Lin-Yau's probability distribution setting with 4 process.
    orf = OllivierRicci(G, alpha=0.5, base=1, exp_power=0, proc=4, verbose="INFO")

    # Do Ricci flow for 2 iterations
    orf.compute_ricci_flow(iterations=2)
    
    orf.set_verbose("ERROR") # mute logs
    orf.compute_ricci_flow(iterations=50)
    G_rf = orf.G.copy()
    
    plt.figure()
    show_results(G_rf)
    plt.savefig("../images/results_rf.png")
    
    plt.figure()
    draw_graph(G_rf)
    plt.savefig("../images/graph_rf.png")
    
    plt.figure()
    draw_graph(my_surgery(G_rf, cut=1.5))
    plt.savefig("../images/graph_rf_surg1.5.png")
    
    plt.figure()
    draw_graph(my_surgery(G_rf, cut=1.0))
    plt.savefig("../images/graph_rf_surg1.0.png")
    
    plt.figure()
    check_accuracy(G_rf)
    plt.savefig("../images/acc_rf.png")
    
    plt.figure()
    draw_graph(my_surgery(G_rf, cut=1.04))
    plt.savefig("../images/graph_rf_surg1.04.png")
    
    orf2 = OllivierRicci(G, alpha=0.5, base=math.e, exp_power=1, proc=8)
    orf2.compute_ricci_flow(iterations=50)
    G_rf2 = orf2.G.copy()
    
    plt.figure()
    show_results(G_rf2)
    plt.savefig("../images/graph_rf_sett2.png")
    
    plt.figure()
    check_accuracy(G_rf2)
    plt.savefig("../images/acc_rf_sett2.png")
    
    plt.figure()
    draw_graph(my_surgery(G_rf2, cut=3)) 
    plt.savefig("../images/graph_rf_surg3.png")
    plt.figure()
    draw_graph(my_surgery(G_rf2, cut=2))
    plt.savefig("../images/graph_rf_surg2.png")
    
    plt.show()
    return    

def check_accuracy(G_origin,weight="weight"):
    """
    Check the accuracy of the Ricci flow edge weight cutoff.
    """
    
    G=G_origin.copy()
    mod, ari = [], []
    
    maxw = max(nx.get_edge_attributes(G, weight).values())
    for cutoff in np.arange(maxw, 0, -0.025):
        trimlist = []
        for n1,n2 in G.edges():
            if G[n1][n2][weight] > cutoff:
                trimlist.append((n1, n2))
        G.remove_edges_from(trimlist)
        cc = list(nx.connected_components(G))
        mod.append(nx.algorithms.community.modularity(G, cc))
        ari.append(ARI(G, cc))
    
    plt.xlim(maxw, 0)
    plt.xlabel("Edge weight cutoff")
    plt.plot(np.arange(maxw,0,-0.025), mod, alpha=0.8)
    plt.plot(np.arange(maxw,0,-0.025), ari, alpha=0.8)
    plt.legend(['Modularity', 'Adjust Rand Index'])
    return

def ARI(G, cc, clustering_label="club"):
    """
    Computer the Adjust Rand Index (clustering accuray) of clustering "cc" with clustering_label as ground truth.
    :param G: A networkx graph
    :param cc: A clustering result as list of connected components list
    :param clustering_label: Node label for clustering groundtruth 
    """
    
    if importlib.util.find_spec("sklearn") is not None:
        from sklearn import preprocessing, metrics
    else:
        print("scikit-learn not installed...")
        return -1
    
    complexlist=nx.get_node_attributes(G, clustering_label)

    le = preprocessing.LabelEncoder()
    y_true=le.fit_transform(list(complexlist.values()))

    predict_dict={}
    for idx, comp in enumerate(cc):
        for c in list(comp):
            predict_dict[c]=idx
    y_pred=[]
    for v in complexlist.keys():
        y_pred.append(predict_dict[v])
    y_pred=np.array(y_pred)
    
    return metrics.adjusted_rand_score(y_true,y_pred)
    
    
def my_surgery(G_origin: nx.Graph(), weight="weight", cut=0):
    """
    A simple surgery function that remove the edges with weight above a threshold
    :param G: A weighted networkx graph
    :param weight: Name of edge weight to cut
    :param cut: Manually assign cut point
    :return: A weighted networkx graph after surgery
    """
    G=G_origin.copy()
    w = nx.get_edge_attributes(G, weight)

    assert cut >=0, "Cut value should be greater than 0."
    if not cut:
        cut = (max(w.values()) - 1.0) * 0.6 + 1.0  # Guess a cut point as default
    
    to_cut = []
    for n1, n2 in G.edges():
        if G[n1][n2][weight] > cut:
            to_cut.append((n1, n2))
    print("*************** Surgery time ****************")
    print("* Cut %d edges." % len(to_cut))
    G.remove_edges_from(to_cut)
    print("* Number of nodes now: %d" % G.number_of_nodes())
    print("* Number of edges now: %d" % G.number_of_edges())
    cc=list(nx.connected_components(G))
    print("* Modularity now: %f " % nx.algorithms.community.modularity(G,cc))
    print("* ARI now: %f " % ARI(G,cc))
    print("*********************************************")
    
    return G    
    
    
def draw_graph(G):
    """
    A helper function to draw a nx graph with community.
    """
    groups =nx.get_node_attributes(G,'club').values()
    color_list = plt.cm.tab10(np.linspace(0, 1, len(groups)))
    color_dict = dict(zip(groups, color_list))
        
    nx.draw_spring(G,seed=0, nodelist=G.nodes(),
                   node_color=[color_dict[x] for x in groups],
                   alpha=0.8)
    return
                   
    
def show_results(G):
    # Print the first five results
    print("Karate Club Graph, first 5 edges: ")
    for n1,n2 in list(G.edges())[:5]:
        print("Ollivier-Ricci curvature of edge (%s,%s) is %f" % (n1 ,n2, G[n1][n2]["ricciCurvature"]))

    # Plot the histogram of Ricci curvatures
    plt.subplot(2, 1, 1)
    ricci_curvtures = nx.get_edge_attributes(G, "ricciCurvature").values()
    plt.hist(ricci_curvtures,bins=20)
    plt.xlabel('Ricci curvature')
    plt.title("Histogram of Ricci Curvatures (Karate Club)")

    # Plot the histogram of edge weights
    plt.subplot(2, 1, 2)
    weights = nx.get_edge_attributes(G, "weight").values()
    plt.hist(weights,bins=20)
    plt.xlabel('Edge weight')
    plt.title("Histogram of Edge weights (Karate Club)")

    plt.tight_layout()
    return    
    
if __name__== "__main__":
    main()  