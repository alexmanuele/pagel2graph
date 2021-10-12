import numpy as np
import pandas as pd

import networkx as nx

def nx_to_dash(G, node):
    nodes = []
    for n in G.nodes:
        if n == node:
            nodes.append({
                        'data': {'id':n, 'label':n, **G.nodes[n]},
                        'classes': 'focal',
            })
        else:
            nodes.append({'data': {'id':n, 'label':n, **G.nodes[n]},
                        'classes':'other',
            })
    edges = []
    for e in G.edges:
        edges.append({'data': {'source': e[0], 'target': e[1], **G.edges[e]}})
    return nodes + edges

def neighborhood(G, node, n):
    path_lengths = nx.single_source_dijkstra_path_length(G, node)
    return [node for node, length in path_lengths.items()
                    if length <= n]

def filter_graph(G, node, d, lr_threshold, p_threshold):
    edges = []
    for u,v,e in G.edges(data=True):
        if e['lr'] >= lr_threshold and e['p'] <= p_threshold:
            edges.append((u,v))
    H=G.edge_subgraph(edges)
    if node in H.nodes:
        return H.subgraph(neighborhood(H, node, d))
    return G.subgraph([node])
