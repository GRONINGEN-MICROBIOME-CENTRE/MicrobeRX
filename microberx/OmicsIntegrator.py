"""
This is a module that provides functions to analyze organims, and enzyme involved in the production of metabolites predicted by MicrobeRX.

The module contains the following functions:

- plot_species_sunburst: The function plot_species_sunburst creates a sunburst plot of the microbial species in the sources list. 

"""

__all__ = ["plot_species_sunburst"]

import plotly.express as px
from plotly.subplots import make_subplots
from distinctipy import distinctipy

import pandas as pd

from .DataFiles import load_microbes_data, load_microbes_reactions


def __check_if_microbes_databases_are_loaded():
    """
    The function check_if_microbes_databases_are_loaded checks if the global variables MICROBES_DATA and MICROBES_REACTIONS are loaded with the data from the microbes databases. If not, it calls the functions load_microbes_data and load_microbes_reactions to load the data and assign it to the global variables.

    Parameters: None

    Returns: None
    """

    global MICROBES_DATA
    MICROBES_DATA = None
    if MICROBES_DATA is not None:
        pass
    else:
        MICROBES_DATA = load_microbes_data()

    global MICROBES_REACTIONS
    MICROBES_REACTIONS = None
    if MICROBES_REACTIONS is not None:
        pass
    else:
        MICROBES_REACTIONS = load_microbes_reactions()


def plot_species_sunburst(sources: list, path: str = "short"):
    """
    The function plot_species_sunburst creates a sunburst plot of the microbial species in the sources list. It uses the global variables MICROBES_DATA and MICROBES_REACTIONS that are loaded by the function check_if_microbes_databases_are_loaded. It also uses the Plotly Express library to create the sunburst plot.

    Parameters
    ----------
        sources: A list of strings that represent the source id of reactions from AGORA2. For example, ['CYSS3r'].
        path: A string that specifies the path of the sunburst plot. It can be either 'full' or 'short'. The default value is 'short'. If 'full', the path is ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']. If 'short', the path is ['Kingdom', 'Phylum', 'Order', 'Genus', 'Species'].
    
    Returns
    -------
    F: A Plotly Figure object that contains the sunburst plot. It has one subplot for each source in the sources list. The subplots are arranged in a grid with three columns and variable rows. The sunburst plot shows the hierarchical distribution of the microbial species by their taxonomic ranks. The color of each segment is determined by the phylum of the species.
    """

    __check_if_microbes_databases_are_loaded()

    selection = MICROBES_REACTIONS[sources].dropna(subset=sources, axis=0, how="all")

    subplots = len(selection.columns)
    cols = 3
    rows = subplots // cols
    if subplots - (rows * cols) > 0:
        rows += 1

    matrix = [dict(row=r + 1, column=c + 1) for r in range(rows) for c in range(cols)]

    specs = [[{"type": "domain"} for c in range(cols)] for r in range(rows)]

    F = make_subplots(
        rows=rows,
        cols=cols,
        horizontal_spacing=0.02,
        vertical_spacing=0.1,
        specs=specs,
        subplot_titles=selection.columns.to_list(),
    )

    phyla = MICROBES_DATA.Phylum.unique()
    colors = [
        distinctipy.get_hex(c)
        for c in distinctipy.get_colors(n_colors=len(phyla), pastel_factor=0.7)
    ]
    color_map = dict(zip(phyla, colors))

    if path == "full":
        path = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    if path == "short":
        path = ["Kingdom", "Phylum", "Order", "Genus", "Species"]

    for index, col in enumerate(selection.columns):
        slide = MICROBES_DATA[
            MICROBES_DATA.microbe_name.isin(selection[[col]].dropna().index)
        ]

        fig1 = px.sunburst(
            slide,
            path=path,
            color="Phylum",
            color_discrete_map=color_map,
            branchvalues="total",
        )

        for d in fig1.data:
            F.add_trace(d, row=matrix[index]["row"], col=matrix[index]["column"])

    F.update_layout(
        grid=dict(
            rows=rows,
            columns=cols,
        ),
        height=400 * rows,
        width=400 * cols,
        title_text=f"Species Sunburst",
    )

    return F
