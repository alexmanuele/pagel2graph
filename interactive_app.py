import numpy as np
import pandas as pd

import networkx as nx

import plotly.graph_objects as go
import plotly.express as px
import plotly.figure_factory as ff

import dendropy

import dash
from dash.dependencies import Output, Input, State
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

# Load extra layouts
cyto.load_extra_layouts()

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

default_stylesheet = [
                        {
                            'selector':'edge',
                            'style': {
                                'width': 'mapData(lr, 50, 200, 0.75, 5)',
                                'opacity': 0.4,
                            },
                        },
                        {'selector': 'node',
                            'style':{
                                #'color': '#317b75',
                                'background-color': '#317b75',
                                'content': 'data(label)',
                            },
                        },
                        {
                            'selector': '.focal',
                            'style':{
                                #'color': '#E65340',
                                'background-color': '#E65340',
                                'content': 'data(label)',
                            },
                        },
                        {'selector': '.other',
                            'style':{
                                #'color': '#317b75',
                                'background-color': '#317b75',
                                'content': 'data(label)',
                            },
                        },
                        ]

def basic_dash_network():
    #Load files for network
    ava_lr = pd.read_table('data/efaecium_profile_LR_rerunNA.csv', sep=',', index_col=0)
    ava_p = pd.read_table('data/efaecium_profile_pval_rerunNA.csv', sep=',', index_col=0)
    ave_lr = pd.read_table('data/pagel_LR_featureVsHabitat.csv', sep=',', index_col=0)
    ave_p = pd.read_table('data/pagel_pvalue_featureVsHabitat.csv', sep=',', index_col=0)

    node_items = [{'label': col, 'value': col} for col in ava_lr.columns]


    G = nx.graphml.read_graphml('data/pagel_results_as_network_updated.graphml')
    node = 'AA893'
    degree = 2
    lr_threshold = 50
    p_threshold = 0.0001
    H = filter_graph(G, node, degree, lr_threshold, p_threshold)
    elements = nx_to_dash(H, 'node')

    #Load files for heatmap thingy
    """
    df = pd.read_table('data/haley_map/Efaecium_ordered_PA_table_Pagel_Aug92021.csv', sep=',')
    df = df.set_index('Isolate')
    #Load the tree, root it, convert to distance matrix
    tree = dendropy.Tree.get(path='data/haley_map/core_gene_alignment.aln.treefile', schema='newick')
    pdm = tree.phylogenetic_distance_matrix()
    table = pdm.as_data_table()
    table.write_csv('dendropy_matrix.csv')
    matr = pd.read_table('dendropy_matrix.csv', sep=',', index_col=0)

    #drop the outgroup since there is no gene info for it.
    dist = matr.drop(list(set(matr.index) - set(df.index)))
    dist = dist.drop(list(set(matr.index) - set(df.index)), axis=1)
    """


    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SOLAR])
    app.layout = dbc.Container(fluid=True, children=[
    #<!-- header -->
        html.H1(children='Visual Analytics Demo'),
        html.Div(children='Explore multiple data sets with Dash.'),
        html.Hr(),

        dbc.Row([
        #<!-- basic plot choices -->
            dbc.Col([
            html.H3(children="Pagel Heatmap"),
            dbc.FormGroup([
                dbc.Label("Choose Metric"),
                dcc.Dropdown(
                    id="dataset-select", value=1,
                    options=[
                        {"label": "Likelihood Ratio", "value": 1},
                        {"label": "P Value", "value": 2},])
                    ]
                ),
            dbc.Button('Update Plot', id='example-button', color='primary', style={'margin-bottom': '1em'}, block=True),
            dbc.Row([
                dbc.Col(dcc.Graph(id='example-graph')),  # Not including fig here because it will be generated with the callback
            ])
        ], lg=3, md=12, className="bg-light text-dark"),

        #<! -- Network plot -- >
        dbc.Col([
            html.H3(children="Network Visualization"),
            dbc.Row([
                dbc.Col(width=6,children=[
                    dbc.FormGroup([
                        dbc.Label("Change network Layout"),
                        dcc.Dropdown(
                            id='network-callbacks-1',
                            value='grid',
                            clearable=False,
                            options=[
                                {'label': name.capitalize(), 'value': name}
                                #for name in ['grid', 'random', 'circle', 'cose', 'concentric', 'breadthfirst']
                                for name in [
                                    'random',
                                    'grid',
                                    'circle',
                                    'concentric',
                                    'breadthfirst',
                                    'cose',
                                    'cose-bilkent',
                                    'cola',
                                    'klay',
                                    'spread',
                                    'euler'
                                ]
                            ], className="bg-light text-dark",
                        ),
                        ]),
                    dbc.FormGroup([
                        dbc.Col([
                            dbc.Label('Select a node of interest.'),
                            dcc.Dropdown(
                                id='node-dropdown',
                                options=node_items,
                                value=node_items[1]['value'],
                                className="bg-light text-dark"),

                        ]),
                        dbc.Col([
                            dbc.Label('Select thresholding values'),
                            dbc.Input(
                                    id='degree',
                                    placeholder='Degree (depth of neighborhood)',
                                    type='number', min=0, step=1, value=2
                                    ),
                            dbc.FormText('Degree (depth of neighborhood)'),
                            dbc.Input(
                                id='lr-threshold',
                                placeholder="Likelihood Ratio lower bound",
                                type='number', min=0, value=50.0
                                ),
                            dbc.FormText('Likelihood Ratio lower bound'),
                            dbc.Input(
                                    id='p-threshold',
                                    placeholder="p-value upper bound",
                                    type='number', min=0, value=0.05),
                            dbc.FormText("p-value upper bound")
                        ]),
                    ]),

                    dbc.Button('Update iPlot', id='interactive-button', color='success', style={'margin-bottom': '1em'}, block=True),
                ]),

                dbc.Col(
                        width=6,
                        children=dbc.Card(
                            [
                                dbc.CardHeader("Network Properties", className="bg-success text-white"),
                                dbc.CardBody(
                                    html.P("Lorem Ipsum and all that.", className='card-text text-dark',
                                    id='node-selected')
                                )
                            ]
                        )
                ),
            ]),
            dbc.Row([
                #dbc.Col(dcc.Graph(id='interactive-graph')),  # Not including fig here because it will be generated with the callback
                dbc.Col(cyto.Cytoscape(
                    id='network-plot',
                    elements=elements,
                    stylesheet=default_stylesheet,
                    style={'width': '100%', 'height': '800px'},
                    layout={
                        'name': 'grid'
                    },
                ),className='bg-white'),

            ]),
        ], lg=9, md=12, className='bg-secondary text-white')
    ]) # </row>

    ],)#</container>


    #plot heatmaps
    @app.callback(
        Output('example-graph', 'figure'),
        [Input('example-button', 'n_clicks'),
        State('dataset-select', 'value'),]
    )
    def plot(click, dataset):
        df_map = {'1': ava_lr,
                  '2': ava_p,}
        plot_map = {'1': px.imshow}
        df = df_map[str(dataset)]
        #plot_f = plot_map[str(plot)]
        plot_f = px.imshow
        return plot_f(df)

    #network layout call back
    @app.callback(
        Output('network-plot', 'layout'),
        Input('network-callbacks-1', 'value')
        )
    def update_layout(layout):
        return {
            'name': layout,
            'animate': True
        }

    @app.callback(
        Output('network-plot', 'elements'),
        Output('node-selected', 'children'),
        [Input('interactive-button', 'n_clicks'),
         State('node-dropdown', 'value'),
         State('degree', 'value'),
         State('lr-threshold', 'value'),
         State('p-threshold', 'value'),]
    )
    def update_elements(click, node, degree, lr_threshold, p_threshold):
        n_nodes = 0
        n_edges = 0
        H = filter_graph(G, node, degree, lr_threshold, p_threshold)
        # Graph basics
        elements = nx_to_dash(H, node)
        n_nodes = len(H.nodes)
        n_edges = len(H.edges)


        #summary = html.P("Focal Node: {0}\nDegree: {1}<br>LR Threshold: {2}<br>P Threshold: {3}<br>Nodes in selection: {4}<br>Edges in selection: {5}".format(node, degree, lr_threshold, p_threshold,n_nodes, n_edges))
        summary = dbc.ListGroup(
            [
                dbc.ListGroupItem("Focal Node: {}".format(node)),
                dbc.ListGroupItem("Degree: {}".format(degree)),
                dbc.ListGroupItem("LR Threshold: {}".format(lr_threshold)),
                dbc.ListGroupItem("P threshold: {}".format(p_threshold)),
                dbc.ListGroupItem("n Nodes: {}".format(n_nodes)),
                dbc.ListGroupItem("n Edges: {}".format(n_edges)),
            ],
        )
        return elements, summary
    @app.callback(Output('network-plot', 'stylesheet'),
                [Input('network-plot', 'tapNode')])
    def highlight_edges(node):
        if not node:
            return default_stylesheet

        stylesheet = [
                            {
                                'selector':'edge',
                                'style': {
                                    'opacity': 0.4,
                                    'width': 'mapData(lr, 50, 200, 0.75, 5)',
                                },
                            },
                            {'selector': 'node',
                                'style':{
                                    #'color': '#317b75',
                                    'background-color': '#317b75',
                                    'content': 'data(label)',
                                    'width': 'mapData(degree, 1, 100, 25, 200)'
                                },
                            },
                            {
                                'selector': '.focal',
                                'style':{
                                    #'color': '#E65340',
                                    'background-color': '#E65340',
                                    'content': 'data(label)',
                                },
                            },
                            {'selector': '.other',
                                'style':{
                                    #'color': '#317b75',
                                    'background-color': '#317b75',
                                    'content': 'data(label)',
                                },
                            },
                            {
                                "selector": 'node[id = "{}"]'.format(node['data']['id']),
                                "style": {
                                    'background-color': '#B10DC9',
                                    "border-color": "purple",
                                    "border-width": 2,
                                    "border-opacity": 1,
                                    "opacity": 1,

                                    "label": "data(label)",
                                    "color": "#B10DC9",
                                    "text-opacity": 1,
                                    "font-size": 12,
                                    'z-index': 9999
                                }
                            }
                        ]
        for edge in node['edgesData']:
            stylesheet.append({
                'selector': 'node[id= "{}"]'.format(edge['target']),
                'style': {
                    'background-color': 'blue',
                    'opacity': 0.9,

                }
            })
            stylesheet.append({
                'selector': 'node[id= "{}"]'.format(edge['source']),
                'style': {
                    'background-color': 'blue',
                    'opacity': 0.9,

                }
            })
            stylesheet.append({
                "selector": 'edge[id= "{}"]'.format(edge['id']),
                "style": {
                    "line-color": 'green',
                    'opacity': 0.9,
                    'z-index': 5000
                }
            })
        return stylesheet
    return app



if __name__ == "__main__":
    app = basic_dash_network()

    app.run_server(debug=True)
