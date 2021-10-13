import itertools as it

import numpy as np
import pandas as pd

import dash
from dash.dependencies import Output, Input, State
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

# Load extra layouts
cyto.load_extra_layouts()

import networkx as nx

import plotly.graph_objects as go
import plotly.express as px

from components import *
from utils import *

ava_lr = pd.read_table('data/efaecium_profile_LR_rerunNA.csv', sep=',', index_col=0)
ava_p = pd.read_table('data/efaecium_profile_pval_rerunNA.csv', sep=',', index_col=0)
ave_lr = pd.read_table('data/pagel_LR_featureVsHabitat.csv', sep=',', index_col=0)
ave_p = pd.read_table('data/pagel_pvalue_featureVsHabitat.csv', sep=',', index_col=0)

node_items = [{'label': col, 'value': col} for col in ava_lr.columns]


G = nx.graphml.read_graphml('data/pagel_results_as_network_updated.graphml')
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
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SOLAR],suppress_callback_exceptions=True)


################################################################################
### Page Layouts                                                             ###
################################################################################

### Entry Point. Also Serves as data dump for sharing between apps  ###
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    #Stores for data persistence.
    dcc.Store(id='graph-store'),

    make_navbar(active=0),
    html.Div(id='page-content'),
])

### Landing Page ###
#Adapted from https://getbootstrap.com/docs/4.0/examples/product/
landing_page_layout = [
    html.Div(className="position-relative overflow-hidden p-3 p-md-5 m-md-3 text-center bg-light",children=[
        html.Div(className="col-md-5 p-lg-5 mx-auto my-5", children=[
            html.H1('ARETE BayesTraits', className="display-4 font-weight-normal"),
            html.P("Interactively explore BayesTraits results.", className='lead font-weight-normal'),
            html.A('Get Started', href='page-1', className='btn btn-outline-secondary'),
        ])
    ]),
    html.Div(className="d-md-flex flex-md-equal w-100 my-md-3 pl-md-3",children=[
        html.Div(className="bg-dark mr-md-3 pt-3 px-3 pt-md-5 px-md-5 text-center text-white overflow-hidden",children=[
            html.Div(className="my-3 py-3", children=[
                html.H2('View LR and p-value heat maps.', className='display-5'),
                html.P(children=["View BayesTraits associations among VF, AMR, GI, and plasmid. Click 'Matrices'.",
                                 ], className='lead'),
            ])
        ]),
        html.Div(className="bg-light mr-md-3 pt-3 px-3 pt-md-5 px-md-5 text-center overflow-hidden",children=[
            html.Div(className="my-3 py-3", children=[
                html.H2("Explore Specific Subnetworks.", className="display-5"),
                html.P('Visualize networks surrounding specific nodes. Choose your node, select filtering parameters, and explore. Click "Network Visualisation."', className='lead'),

            ])
        ])
    ]),
    html.Div(className="d-md-flex flex-md-equal w-100 my-md-3 pl-md-3",children=[
        html.Div(className="bg-light mr-md-3 pt-3 px-3 pt-md-5 px-md-5 text-center overflow-hidden",children=[
            html.Div(className="my-3 py-3", children=[
                html.H2('View Macro Level Statistics.', className='display-5'),
                html.P("View at a macro level how LR and p value choices influence network properties. Click 'Network Statistics'.", className='lead'),
            ])
        ]),
        html.Div(className="bg-primary mr-md-3 pt-3 px-3 pt-md-5 px-md-5 text-center text-white overflow-hidden",children=[
            html.Div(className="my-3 py-3", children=[
                html.H2("(Not yet) View genome presence/absence cladogram.", className="display-5"),
                html.P('See the presence and absence of genes etc. in a per genome basis. (Not yet implemented).', className='lead'),

            ])
        ])
    ]),
]

### Heat Map Viewer Layout ###
page1_layout = dbc.Container(fluid=True, children=[
    dbc.Row(id='heatmap-display',children=[
        dbc.Col([
            dcc.Loading(dcc.Graph(id='heatmap-graph'),),
            dbc.Row([
                dbc.Col([
                    dbc.FormGroup([
                        dbc.Label("Choose Metric"),
                        dcc.Dropdown(
                            id="dataset-select", value=1,
                            options=[
                                {"label": "Likelihood Ratio", "value": 1},
                                {"label": "P Value", "value": 2},])
                            ]
                        ),
                ],className='pl-5 pr-5'),
            ],),
        ]),
    ]),
])


### Network summary layout ###
page2_layout = dbc.Container(fluid=True, children=[
    dbc.Row([
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
                    elements=[],
                    stylesheet=default_stylesheet,
                    style={'width': '100%', 'height': '800px'},
                    layout={
                        'name': 'grid'
                    },
                ),className='bg-white'),

            ]),
        ], className='bg-secondary text-white')
    ]),
])

page3_layout = dbc.Container(fluid=True, children=[
    dbc.Row(id='historgram-display',children=[
        dbc.Col([
            dcc.Loading(dcc.Graph(id='histogram-graph'),),
            dbc.Row([
                dbc.Col([
                    dbc.FormGroup([
                        dbc.Label("Choose Facet Metric"),
                        dcc.Dropdown(
                            id="histogram-metric-select", value=1,
                            options=[
                                {"label": "Likelihood Ratio", "value": 1},
                                {"label": "P Value", "value": 2},]),
                        dbc.Label("Choose y Value"),
                        dcc.Dropdown(
                            id="histogram-y-select", value=1,
                            options=[
                                {"label": "Focal Node Degree", "value": 1},
                                {"label": "Graph n Nodes", "value": 2},
                                {'label': "Graph n Edges", 'value': 3},]),
                        dbc.Button('Re-calculate Plot', id='histogram-button', color='primary', style={'margin-bottom': '1em'}, block=True),
                    ]),
                ],className='pl-5 pr-5'),
            ],),
        ]),
    ]),
])



################################################################################
### Heatmap Callbacks                                                        ###
################################################################################
@app.callback(
    Output('heatmap-graph', 'figure'),
    [Input('dataset-select', 'value'),]
)
def plot(dataset):
    df_map = {'1': (ava_lr, ave_lr),
              '2': (ava_p, ave_p),}

    colorbar_map = {'1': dict(
                    tick0=0,
                    dtick=50,
                    title='LR',
                    tickvals=[0, 50, 100,150,200],
                    ticktext=['< 0', '50', '100', '150', '200+'],
                    tickmode='array',
                    ),
                    '2': dict(tick0=0,
                              dtick=0.05,
                              tickvals=[0.0, 0.05, 0.15, 1.0],
                              tickmode='array',
                              title='p-value',),
                    }

    colorscale_map = {'1': [(0.00, px.colors.sequential.ice[0]),   (0.03846, px.colors.sequential.ice[0]),
                            (0.03846, px.colors.sequential.ice[2]), (0.23076923076923078, px.colors.sequential.ice[2]),
                            (0.23076923076923078, px.colors.sequential.ice[4]),  (0.42307692307692313, px.colors.sequential.ice[4]),
                            (0.42307692307692313, px.colors.sequential.ice[6]), (0.6153846153846154, px.colors.sequential.ice[6]),
                            (0.6153846153846154, px.colors.sequential.ice[8]), (0.8076923076923077, px.colors.sequential.ice[8]),
                            (0.8076923076923077, px.colors.sequential.ice[10]), (1.0, px.colors.sequential.ice[10])],
                      '2': [(0.00, px.colors.sequential.ice[0]),   (0.05, px.colors.sequential.ice[0]),
                               (0.05, px.colors.sequential.ice[5]), (0.15, px.colors.sequential.ice[3]),
                               (0.15, px.colors.sequential.ice[-1]),  (1.00, px.colors.sequential.ice[-1])]}

    z_bounds_map = {'1': (-10, 250),
                    '2': (0, 1.0)}

    ava, ave = df_map[str(dataset)]
    zmin, zmax = z_bounds_map[str(dataset)]
    colorbar = colorbar_map[str(dataset)]
    colorscale = colorscale_map[str(dataset)]

    fig = go.Figure()

    ave_hm = go.Heatmap(x=ave.columns,
                        y=ave.index,
                        z=ave,
                        zmin=zmin,
                        zmax=zmax,
                        colorscale=colorscale,
                        showscale=False,
                       )

    ava_hm = go.Heatmap(x=ava.columns,
                        y=ava.index,
                        z=ava,
                        zmin=zmin,
                        zmax=zmax,
                        colorscale=colorscale,
                        colorbar=colorbar,
                       )

    f1 = go.Figure(ave_hm)
    for data in f1.data:
        fig.add_trace(data)

    f2 = go.Figure(ava_hm)
    for i in range(len(f2['data'])):
        f2['data'][i]['xaxis'] = 'x2'

    for data in f2.data:
        fig.add_trace(data)

    fig.update_layout({'height':800})
    fig.update_layout(xaxis={'domain': [.0, .20],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      #'ticks':"",
                                      #'showticklabels': False
                            })
    # Edit xaxis2
    fig.update_layout(xaxis2={'domain': [0.25, 1.0],
                                       'mirror': False,
                                       'showgrid': False,
                                       'showline': False,
                                       'zeroline': False,
                                       #'showticklabels': False,
                                       #'ticks':""
                             })

    return fig

################################################################################
### Network Visualization Callbacks                                          ###
################################################################################

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

################################################################################
### Network Visualization Callbacks                                          ###
################################################################################
@app.callback(
    Output('histogram-graph', 'figure'),
    [Input('histogram-button', 'n_clicks'),
    State('histogram-metric-select', 'value'),
    State('histogram-y-select', 'value')]
)
def show_histogram(click, metric_sel, y_sel):
    metric_map = {
                    '1': ('lr', 'p'),
                    '2': ('p', 'lr')
                }
    y_map = {
                '1': 'node_degree',
                '2': 'n_nodes',
                '3': 'n_edges'
            }
    dynamic_metric, static_metric = metric_map[str(metric_sel)]
    y = y_map[str(y_sel)]


    lte = lambda x,y: x<=y
    gte = lambda x,y: x>=y

    if dynamic_metric == 'lr':
        search = [25, 50, 100, 150]
        dfun = gte
        sfun = lte

    else:
        search = [0.05, 1e-5, 1e-9, 1e-12]
        dfun = lte
        sfun = gte

    if static_metric == 'p':
        static_threshold = 0.05
    else:
        static_threshold = 50

    records = []
    node_list = []
    for node in G.nodes:
        node_list.append((node, {**G.nodes[node]}))

    for dynamic_threshold in search:
        F= nx.Graph()
        F.add_nodes_from(node_list)

        edges = []
        for u,v,e in G.edges(data=True):
            if dfun(e[dynamic_metric], dynamic_threshold) and sfun(e[static_metric], static_threshold):
            #if e['lr'] >= lr_threshold and e['p'] <= p_threshold:
                edges.append((u,v, e))

        F.add_edges_from(edges)

        graph_degree = F.degree()
        for i, node in enumerate(F.nodes):
            Sub = F.subgraph(neighborhood(F, node, 2))
            n_nodes = len(Sub.nodes)
            n_edges = len(Sub.edges)

            records.append({'node': node,
                            'node_degree': graph_degree[node],
                            'n_nodes': n_nodes,
                            'n_edges': n_edges,
                            dynamic_metric: dynamic_threshold,
                            static_metric: static_threshold,})

    rdf = pd.DataFrame.from_records(records)
    plot = px.histogram(rdf, x='node', y=y, facet_col=dynamic_metric)
    plot.update_layout({'height':800})
    return plot


################################################################################
### Page Navigation callbacks                                                ###
################################################################################
@app.callback(
    [Output('page-content', 'children'),
    Output('page-1-nav', 'className'),
    Output('page-2-nav', 'className'),
    Output('page-3-nav', 'className'),],
    [Input('url', 'pathname'),]
)
def display_page(pathname):
    if pathname == '/page-1':
        return page1_layout, 'active', '', '',
    elif pathname == '/page-2':
        return page2_layout, '', 'active', '',

    elif pathname == '/page-3':
        return page3_layout, '', '', 'active',

    else:
        return landing_page_layout, '', '', '',


if __name__ == '__main__':
    app.run_server(debug=True)
