# Some components are reused in each app. Put here for easier code readability
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

def make_navbar(active=0):
    classnames = ['', '', '']
    classnames[active] = "active"

    navbar = dbc.NavbarSimple(
        children=[
            dbc.NavItem(dbc.NavLink("Matrices", href="/page-1"),id='page-1-nav' ,className=classnames[0]),
            dbc.NavItem(dbc.NavLink("Network Visualization", href="page-2"),id='page-2-nav', className=classnames[1]),
            dbc.NavItem(dbc.NavLink("Network Statistics", href="page-3"),id='page-3-nav', className=classnames[2]),
        ],
        brand="ARETE BayesTraits",
        brand_href="/",
        color="primary",
        dark=True,
    )
    return navbar
