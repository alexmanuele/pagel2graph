import numpy as np
import pandas as pd

import networkx as nx

import plotly.graph_objects as go
import plotly.express as px

import dash
from dash.dependencies import Output, Input, State
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto


def nx_to_dash(G):
    nodes = []
    for n in G.nodes:
        nodes.append({'data': {'id':n, 'label':n, **G.nodes[n]}})
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

def basic_dash_network():
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
    # Graph basics
    elements = nx_to_dash(H)


    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
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
                dcc.Dropdown(id="dataset-select", value=1,
                             options=[{"label": "Likelihood Ratio", "value": 1},
                                      {"label": "P Value", "value": 2},])

                ]),

            dbc.Button('Update Plot', id='example-button', color='primary', style={'margin-bottom': '1em'}, block=True),
            dbc.Row([
                dbc.Col(dcc.Graph(id='example-graph')),  # Not including fig here because it will be generated with the callback
            ])
        ], lg=6, md=12, className="bg-light text-dark"),

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
                                for name in ['grid', 'random', 'circle', 'cose', 'concentric']
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
                            dbc.Input(id='degree',
                                     placeholder='Degree (depth of neighborhood)',
                                     type='number', min=0, step=1, value=2),
                            dbc.FormText('Degree (depth of neighborhood)'),
                            dbc.Input(id='lr-threshold',
                                      placeholder="Likelihood Ratio lower bound",
                                      type='number', min=0, value=50.0),
                            dbc.FormText('Likelihood Ratio lower bound'),
                            dbc.Input(id='p-threshold',
                                      placeholder="p-value upper bound",
                                      type='number', min=0, value=0.05),
                            dbc.FormText("p-value upper bound")
                        ]),
                    ]),

                    dbc.Button('Update iPlot', id='interactive-button', color='success', style={'margin-bottom': '1em'}, block=True),
                ]),

                dbc.Col(width=6,
                  children=dbc.Card(
                  [
                    dbc.CardHeader("Network Properties", className="bg-success text-white"),
                    dbc.CardBody(
                        html.P("Lorem Ipsum and all that.", className='card-text text-dark',
                        id='node-selected')
                    )

                  ])
                ),
            ]),
            dbc.Row([
                #dbc.Col(dcc.Graph(id='interactive-graph')),  # Not including fig here because it will be generated with the callback
                dbc.Col(cyto.Cytoscape(
                    id='network-plot',
                    elements=elements,
                    style={'width': '100%', 'height': '400px'},
                    layout={
                        'name': 'grid'
                    }
                )),

            ]),
        ], lg=6, md=12, className='bg-secondary text-white')
    ]) # </row>

    ],)#</container>


    #Make the plot for the 3 datasets (plot 1, i suppose)
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
        elements = nx_to_dash(H)
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
    return app


if __name__ == "__main__":
    app = basic_dash_network()

    app.run_server(debug=True)
