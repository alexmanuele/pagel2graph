import networkx as nx
import argparse

argparser = argparse.ArgumentParser(description='Filter GraphML file to explore relationships.')
requiredNamed = argparser.add_argument_group('required named arguments')
requiredNamed.add_argument('-i', help='Input GraphML file.', required=True)
requiredNamed.add_argument('-n', help='Node of interest. Must be an exact match with a node in the graph.', required=True)
requiredNamed.add_argument('-d', help='Degree of neighborhood from node of interest to include.', required=True, type=int)
requiredNamed.add_argument('-lr', help='Likelihood ratio threshold. Edges below this value will be excluded.', required=True, type=float)
requiredNamed.add_argument('-p', help='P-value ratio threshold. Edges above this value will be excluded.', required=True, type=float)
requiredNamed.add_argument('-o', help='Path for output file.', required=True)

#Function uses dijkstra to calculate paths through the graph.
#Given a node and a degree, returns the nodes within degree n
def neighborhood(G, node, n):
    path_lengths = nx.single_source_dijkstra_path_length(G, node)
    return [node for node, length in path_lengths.items()
                    if length <= n]

if __name__=='__main__':
    args = argparser.parse_args()

    inpath = args.i
    outpath = args.o
    node = args.n
    degree = args.d
    lr_threshold = args.lr
    p_threshold = args.p

    G = nx.graphml.read_graphml(inpath)

    try:
        assert node in G.nodes
    except:
        print("Node {} was not found in the graph. Please double check spelling of the node and file path.")
        exit()

    edges = []
    for u,v,e in G.edges(data=True):
        if e['lr'] >= lr_threshold and e['p'] <= p_threshold:
            edges.append((u,v))

    H = G.edge_subgraph(edges)

    try:
        selected = neighborhood(H, node, degree)
    except:
        print("The node was not found in the filtered graph.")
        print("Try specifying a different node or different thresholds.")
        exit()

    nx.readwrite.graphml.write_graphml(H.subgraph(selected), outpath)
