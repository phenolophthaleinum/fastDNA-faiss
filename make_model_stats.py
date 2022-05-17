import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output, dash_table
import plotly.express as px


# df_f = pd.read_pickle("df_stats.pkl")
# fig = px.bar(df_f, x='names', y='genome counts', color="taxlevel", barmode="group")
# fig.show()

# df_g = pd.read_pickle("df_general.pkl")
# df_val = df_g == ''
# df_val_final = df_val.iloc[:, [1, 2, 3, 4, 5]]
# heat_cells = df_val_final.to_numpy().astype(int)
# heat_cells_annot = df_g.iloc[:, [1, 2, 3, 4, 5]].to_numpy()
# x_heat = list(df_g.iloc[:, [1, 2, 3, 4, 5]])
# y_heat = df_g['organism'].tolist()
# my_colorsc = [[0, 'rgb(0,255,0)'],
#               [1, 'rgb(255,0,0']]
# fig = px.imshow(heat_cells, x=x_heat, y=y_heat,
#                 color_continuous_scale=[(0.00, "rgb(64, 237, 90)"), (1.00, "rgb(237, 64, 64)")], aspect="auto")
# fig.update_traces(text=heat_cells_annot, hovertemplate="%{text}")
# fig.update_xaxes(side="top")
# fig.show()

app = Dash(
    external_stylesheets=[dbc.themes.SANDSTONE],
)
df_final = pd.read_pickle("df_stats.pkl")
df_g = pd.read_pickle("df_general.pkl")
# app.config.suppress_callback_exceptions = True

# controls = dbc.Card(
#     [
#         html.Div(
#             [
#                 dbc.Label("Measurement (Y-axis)"),
#                 dcc.Dropdown(
#                     id="dropdown",
#                     options=['genome counts', 'genome length avg', 'genome length std'],
#                     value="genome counts",
#                     clearable=False,
#                 ),
#
#             ]
#         ),
#         html.Div(
#             [
#                 dbc.Label("Plot height"),
#                 dcc.Slider(
#                     id="h-slider", min=400, max=1200, step=200, value=400,
#                     marks={x: str(x) for x in range(400, 1200, 200)},
#                 ),
#             ]
#         ),
#     ],
#     body=True,
# )
# df_final = pd.read_pickle("df_stats.pkl")
# app.layout = html.Div([
#     html.H4('make-model stats'),
#     dcc.Dropdown(
#         id="dropdown",
#         options=['genome counts', 'genome length avg', 'genome length std'],
#         value="genome counts",
#         clearable=False,
#     ),
#     dcc.Graph(id="graph"),
# ])

# app.layout = dbc.Container(
#     [
#         html.H1("make-model stats"),
#         html.Hr(),
#         dcc.Tabs([
#             dcc.Tab(label='General stats', children=[
#                 dbc.Row(
#                     [
#                         controls,
#                         dcc.Graph(id="graph"),
#                     ],
#                     align="center",
#                 ),
#             ]),
#             dcc.Tab(label='Missing names', children=[
#                 dcc.Graph(id="heat-graph"),
#             ]),
#         ]),
#     ],
#     fluid=True,
# )
app.layout = dbc.Container(
    [
        dcc.Store(id="store"),
        html.H1("make-model stats"),
        html.Hr(),
        dbc.Card(
            [
                html.Div(
                    [
                        dbc.Label("Measurement (Y-axis)"),
                        dcc.Dropdown(
                            id="dropdown",
                            options=['genome counts', 'genome length avg', 'genome length std'],
                            value="genome counts",
                            clearable=False,
                        ),

                    ]
                ),
                html.Div(
                    [
                        dbc.Label("Plot height"),
                        dcc.Slider(
                            id="h-slider", min=400, max=1200, step=200, value=400,
                            marks={x: str(x) for x in range(400, 1200, 200)},
                        ),
                    ]
                ),
            ],
            body=True,
        ),
        dbc.Tabs(
            [
                dbc.Tab(label="General stats", tab_id="general"),
                dbc.Tab(label="Missing names", tab_id="heat"),
            ],
            id="tabs",
            active_tab="scatter",
        ),
        html.Div(id="tab-content", className="p-4"),
    ],
    fluid=True,
)


@app.callback(
    Output("tab-content", "children"),
    [Input("tabs", "active_tab"), Input("store", "data")],
)
def render_tab_content(active_tab, data):
    """
    This callback takes the 'active_tab' property as input, as well as the
    stored graphs, and renders the tab content depending on what the value of
    'active_tab' is.
    """
    if active_tab and data is not None:
        if active_tab == "general":
            return dbc.Row(
                [
                    dash_table.DataTable(
                        id='datatable-interactivity',
                        columns=[
                            {"name": i, "id": i, "deletable": False, "selectable": True} for i in df_final.columns
                        ],
                        data=df_final.to_dict('records'),
                        editable=True,
                        filter_action="native",
                        sort_action="native",
                        sort_mode="multi",
                        column_selectable="single",
                        row_selectable="multi",
                        row_deletable=False,
                        selected_columns=[],
                        selected_rows=[],
                        page_action="native",
                        page_current=0,
                        page_size=10,
                    ),
                    dcc.Graph(figure=data['general'])
                ],
                align="center",
            ),
        elif active_tab == "heat":
            return dcc.Graph(figure=data['heat'])
    return "No tab selected"


# @app.callback(
#     Output("graph", "figure"),
#     [
#         Input("dropdown", "value"),
#         Input("h-slider", "value")
#     ])
# def update_bar_chart(measure, height):
#     df_final = pd.read_pickle("df_stats.pkl")  # replace with your own data source
#     # mask = df["day"] == measure
#     x_data = [df_final["taxlevel"].tolist(), df_final["names"].tolist()]
#     fig = go.Figure()
#     fig.add_bar(x=x_data, y=df_final[measure])
#     # fig = px.bar(df_final, x=x_data, y=measure,
#     #              color="taxlevel", barmode="group")
#     fig.update_layout(height=int(height), barmode='group')
#     return fig
@app.callback(
    Output("store", "data"),
    [
        Input("dropdown", "value"),
        Input("h-slider", "value")
    ])
def generate_graphs(measure, height):
    """
    This callback generates three simple graphs from random data.
    """
    df_final = pd.read_pickle("df_stats.pkl")  # replace with your own data source
    # mask = df["day"] == measure
    x_data = [df_final["taxlevel"].tolist(), df_final["names"].tolist()]
    fig_bar = go.Figure()
    fig_bar.add_bar(x=x_data, y=df_final[measure])
    # fig = px.bar(df_final, x=x_data, y=measure,
    #              color="taxlevel", barmode="group")
    fig_bar.update_layout(height=int(height), barmode='group')

    df_g = pd.read_pickle("df_general.pkl")
    df_val = df_g == ''
    df_val_final = df_val.iloc[:, [1, 2, 3, 4, 5]]
    heat_cells = df_val_final.to_numpy().astype(int)
    heat_cells_annot = df_g.iloc[:, [1, 2, 3, 4, 5]].to_numpy()
    x_heat = list(df_g.iloc[:, [1, 2, 3, 4, 5]])
    y_heat = df_g['organism'].tolist()
    my_colorsc = [[0, 'rgb(0,255,0)'],
                  [1, 'rgb(255,0,0']]
    fig_heat = px.imshow(heat_cells, x=x_heat, y=y_heat,
                         color_continuous_scale=[(0.00, "rgb(64, 237, 90)"), (1.00, "rgb(237, 64, 64)")], aspect="auto")
    fig_heat.update_traces(text=heat_cells_annot, hovertemplate="%{text}")
    fig_heat.update_xaxes(side="top")
    fig_heat.update_layout(height=int(height))
    # save figures in a dictionary for sending to the dcc.Store
    return {"general": fig_bar, "heat": fig_heat}


app.run_server(debug=True)
# fig = px.bar(df_final, x='names', y='genome counts', color="taxlevel", barmode="group")
# fig.show()
