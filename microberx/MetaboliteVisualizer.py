"""
This is a module that provides functions to visualize predicted metabolites from MicrobeRX.

The module contains the following functions:

- plot_molecular_descriptors: Generated and interactive of the molecular descriptors of a given data frame using polar coordinates.

- plot_isotopic_masses: Generated and interactive plot of the isotopic mass distribution of a given data frame using plotly.

- plot_confidence_scores: Creates a 3D scatter plot of the data frame with the x, y, and z axes representing the similarity of substrates, products, and reacting atoms efficiency respectively.

- plot_metabolic_accesibility: Creates a 2D image of a molecule with the atoms colored according to their metabolic accessibility.

- plot_relationships: Creates a Sankey diagram to visualize the evidences of metabolite annotations in a data frame.

- display_molecules: Displays a grid of molecules from a data frame, using different colors to indicate the values of a specified column.

"""

__all__ = [
    "plot_molecular_descriptors",
    "plot_isotopic_masses",
    "plot_confidence_scores",
    "plot_metabolic_accesibility",
    "display_molecules",
    "plot_relationships",
]

import io, copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

import plotly.express as px
import plotly.graph_objects as go

from distinctipy import distinctipy

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

from PIL import Image

import mols2grid

import warnings

warnings.filterwarnings("ignore")


def plot_molecular_descriptors(data_frame: pd.DataFrame, names_col: str):
    """
    Plots the molecular descriptors of a given data frame using polar coordinates.

    Parameters
    ----------
    data_frame : pd.DataFrame
        A pandas data frame that contains the molecular descriptors as columns and the compound names as rows.
    names_col : str
        A string that specifies the name of the column that contains the compound names.

    Returns
    -------
    Figure : plotly.graph_objects.Figure
        A plotly figure object that shows the polar plot of the molecular descriptors. The plot has the following features:
            - The radial axis represents the normalized value of each molecular descriptor, ranging from 0 to 1.
            - The angular axis represents the different molecular descriptors, such as MolWt, LogP, NumHAcceptors, etc.
            - Each compound is plotted as a radial line with a distinct color and a marker at each descriptor value.
            - The upper and lower limits of the Lipinski's rule of five are plotted as shaded regions in orange and yellow, respectively. The rule of five states that most drug-like molecules have molecular weight less than 500, LogP less than 5, number of hydrogen bond acceptors less than 10, and number of hydrogen bond donors less than 5.
            - A legend is displayed on the right side of the plot, showing the name and color of each compound.
    """

    data = pd.DataFrame(index=list(data_frame[names_col]))
    data["MolWt"] = [i / 500 for i in data_frame["MolWt"]]
    data["LogP"] = [i / 5 for i in data_frame["LogP"]]
    data["HBA"] = [i / 10 for i in data_frame["NumHAcceptors"]]
    data["HBD"] = [i / 5 for i in data_frame["NumHDonors"]]
    data["RotB"] = [i / 10 for i in data_frame["NumRotatableBonds"]]
    data["TPSA"] = [i / 140 for i in data_frame["TPSA"]]

    Ro5_up = [1, 1, 1, 1, 1, 1]
    Ro5_low = [0.5, 0.1, 0.1, 0.25, 0.1, 0.5]
    # data=data.reindex(natsorted(data.index))
    categories = list(data.columns)

    Figure = go.Figure()

    fig1 = px.line_polar(
        r=Ro5_up,
        theta=categories,
        line_close=True,
        color_discrete_sequence=["red"],
    )
    fig1.update_traces(
        fill="toself", fillcolor="orange", opacity=0.4, name="Upper limit"
    )
    fig2 = px.line_polar(
        r=Ro5_low, theta=categories, line_close=True, color_discrete_sequence=["red"]
    )
    fig2.update_traces(
        fill="toself", fillcolor="yellow", opacity=0.4, name="Lower limit"
    )

    fig1.data[0].showlegend = True
    fig2.data[0].showlegend = True

    Figure.add_trace(fig1.data[0])
    Figure.add_trace(fig2.data[0])

    color = [
        distinctipy.get_hex(c)
        for c in distinctipy.get_colors(n_colors=len(data.index), pastel_factor=0.7)
    ]

    for i, name in enumerate(data.index):
        f = px.line_polar(
            r=list(data.loc[name]),
            theta=categories,
            line_close=True,
            color_discrete_sequence=[color[i]],
        )
        f.update_traces(name=name)
        f.data[0].showlegend = True
        Figure.add_trace(f.data[0])

    Figure.update_layout(
        width=800,
        height=600,
        polar=dict(radialaxis=dict(visible=True)),
        showlegend=True,
        legend=dict(x=1.2, y=0.95),
        legend_orientation="h",
    )

    return Figure


def plot_isotopic_masses(data_frame: pd.DataFrame, names_col: str, mass_distribution_col: str):
    """
    Plots the isotopic mass distribution of a given data frame using plotly.

    Parameters
    ----------
    data_frame : pd.DataFrame
        A pandas data frame that contains the isotopic mass distribution as a column of strings, where each string has the format 'mass:probability;mass:probability;...'
    names_col : str
        A string that specifies the name of the column that contains the compound names.
    mass_distribution_col : str
        A string that specifies the name of the column that contains the isotopic mass distribution.

    Returns
    -------
    Figure : plotly.graph_objects.Figure
        A plotly figure object that shows the bar plot of the isotopic mass distribution for each compound. The plot has the following features:
            - The x-axis represents the mass values of the isotopes, rounded to four decimal places.
            - The y-axis represents the probability values of the isotopes, multiplied by 100 and rounded to four decimal places.
            - Each compound is plotted as a group of bars with a distinct color and a label at the top of each bar.
            - A legend is displayed on the right side of the plot, showing the name and color of each compound.
    """

    data_frame = copy.deepcopy(data_frame)
    masses = data_frame[[names_col, mass_distribution_col]]
    masses.dropna(subset=[mass_distribution_col], inplace=True)
    masses.mass_distribution = masses[mass_distribution_col].apply(
        lambda x: {
            float(value.split(":")[0]): float(value.split(":")[1])
            for value in x.split(";")
        }
    )
    masses.reset_index(drop=True, inplace=True)

    Figure = go.Figure()

    color = [
        distinctipy.get_hex(c)
        for c in distinctipy.get_colors(n_colors=len(masses.index), pastel_factor=0.7)
    ]

    for index in masses.index:
        f = px.bar(
            x=tuple(masses[mass_distribution_col][index].keys()),
            y=tuple(masses[mass_distribution_col][index].values()),
            color_discrete_sequence=[color[index]],
        )
        f.update_traces(name=masses[names_col][index])
        f.data[0].showlegend = True
        Figure.add_trace(f.data[0])

    Figure.update_layout(
        width=800,
        height=600,
        showlegend=True,
        legend=dict(x=1.2, y=0.95),
        legend_orientation="h",
        plot_bgcolor="rgba(0, 0, 0, 0)",
        font=dict(size=12),
        xaxis_title="Atomic Mass (u)",
        yaxis_title="Relative abundance (%)",
        title="Isotopic distribution",
    )

    return Figure


def plot_confidence_scores(data_frame: pd.DataFrame,x: str = "similarity_substrates",y: str = "similarity_products",z: str = "reacting_atoms_efficiency",cmap: str = "RdYlGn"):
    """
    Creates a 3D scatter plot of the data frame with the x, y, and z axes representing the similarity of substrates, products, and reacting atoms efficiency respectively.

    Parameters
    ----------
    data_frame : pd.DataFrame
        The data frame containing the columns ‘similarity_substrates’, ‘similarity_products’, ‘reacting_atoms_efficiency’, ‘confidence_score’, and ‘metabolite_id’.
    x : str, optional
        The name of the column to use as the x-axis. Defaults to ‘similarity_substrates’.
    y : str, optional
        The name of the column to use as the y-axis. Defaults to ‘similarity_products’.
    z : str, optional
        The name of the column to use as the z-axis. Defaults to ‘reacting_atoms_efficiency’.
    cmap : str, optional
        The name of the color map to use for the color scale. Defaults to ‘RdYlGn’.

    Returns
    -------
    plotly.Figure
        The 3D scatter plot figure. The figure has the following features:
            - The x-axis represents the similarity of substrates, ranging from 0 to 1.
            - The y-axis represents the similarity of products, ranging from 0 to 1.
            - The z-axis represents the reacting atoms efficiency, ranging from 0 to 1.
            - The color of each point indicates the confidence score of the corresponding metabolite id, ranging from 0 to 1. A color bar is displayed on the right side of the plot.
            - The hover text of each point shows the metabolite id and the values of x, y, z, and color.
            - The title of the plot shows the names of the columns used for x, y, z, and color.
    """
    fig = px.scatter_3d(
        data_frame,
        x="similarity_substrates",
        y="similarity_products",
        z="reacting_atoms_efficiency",
        color="confidence_score",
        range_color=[0, 3],
        color_continuous_scale=cmap,
        width=800,
        height=600,
        opacity=0.0,
        hover_name="metabolite_id",
    )

    fig.update_traces(
        marker_size=5,
        showlegend=False,
    )

    fig.update_layout(
        paper_bgcolor="white",
        plot_bgcolor="white",
        scene=dict(
            xaxis=dict(
                title="Similarity substrates",
                showgrid=True,
                showline=True,
                linecolor="black",
                backgroundcolor="rgb(200, 200, 230)",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="black",
                range=(0, 1),
                dtick=0.25,
                nticks=4,
            ),
            yaxis=dict(
                title="Similarity products",
                showgrid=True,
                showline=True,
                linecolor="black",
                backgroundcolor="rgb(230, 200,230)",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="black",
                range=(0, 1),
                dtick=0.25,
                nticks=4,
            ),
            zaxis=dict(
                title="Atoms efficiency",
                showgrid=True,
                showline=True,
                linecolor="black",
                backgroundcolor="rgb(230, 230,200)",
                gridcolor="white",
                showbackground=True,
                zerolinecolor="black",
                range=(0, 1),
                dtick=0.25,
                nticks=4,
            ),
        ),
        coloraxis_colorbar_title_text="Confidence score",
        coloraxis_colorbar_len=0.5,
        coloraxis_colorbar_dtick=0.5,
    )
    return fig


def plot_metabolic_accesibility(data_frame: pd.DataFrame,molecule: Chem.Mol,atom_map_col: str = "reacting_atoms_in_query",mol_name: str = "Query",alpha: float = 0.5,cmap: str = "RdYlGn_r"):
    """
    Creates a 2D image of a molecule with the atoms colored according to their metabolic accessibility.

    Parameters
    ----------
    data_frame : pd.DataFrame
        The data frame containing the column with the atom map information.
    molecule : Chem.Mol
        The molecule object to be drawn.
    atom_map_col : str, optional
        The name of the column with the atom map information. The column should contain lists of integers representing the atom indices. Defaults to 'reacting_atoms_in_query'.
    mol_name : str, optional
        The name of the molecule to be displayed on the image. Defaults to 'Query'.
    alpha : float, optional
        The transparency level of the atom colors, ranging from 0 to 1. Defaults to 0.5.
    cmap : str, optional
        The name of the color map to use for the color scale. Defaults to 'RdYlGn_r'.

    Returns
    -------
    matplotlib.Figure
        The 2D image figure. The figure has the following features:
            - The molecule is drawn in a 2D projection with the atom symbols and bond types shown.
            - The atoms are colored according to their metabolic accessibility, which is calculated as the frequency of the atom in the atom map column of the data frame. The color scale ranges from red (low accessibility) to green (high accessibility).
            - A color bar is displayed on the right side of the image, showing the values of the metabolic accessibility.
            - The name of the molecule is displayed on the top left corner of the image.
    """

    reacting_atoms = [a for l in data_frame[atom_map_col] for a in l]
    pyplot_cmap = plt.get_cmap(cmap)
    atoms_counts = {i: reacting_atoms.count(i) for i in reacting_atoms}

    def _scale_value(val, values):
        """
        Scales the value between 0 and 1.
        """
        return (val - min(values)) / (max(values) - min(values))

    color_atom_map = {
        k: [pyplot_cmap(_scale_value(v, atoms_counts.values()), alpha=alpha)]
        for k, v in atoms_counts.items()
    }
    d2d = rdMolDraw2D.MolDraw2DCairo(1600, 400)
    d2d.ClearDrawing()
    dos = d2d.drawOptions()
    dos.bondLineWidth = 8
    dos.centreMoleculesBeforeDrawing = True
    dos.highlightRadius = 0.5
    dos.useBWAtomPalette()
    d2d.DrawMoleculeWithHighlights(
        molecule,
        legend="",
        highlight_atom_map=color_atom_map,
        highlight_bond_map={},
        highlight_radii={},
        highlight_linewidth_multipliers={},
    )

    d2d.FinishDrawing()
    img = Image.open(io.BytesIO(d2d.GetDrawingText()))

    fig, axes = plt.subplots(figsize=(8, 15))
    ax = plt.imshow(img)
    ax.axes.set_title(mol_name, fontsize=18, fontweight="bold")
    ax.axes.set_xticklabels([])
    ax.axes.set_yticklabels([])
    ax.axes.tick_params(axis="both", which="both", length=0)
    for spine in ["bottom", "top", "right", "left"]:
        ax.axes.spines[spine].set_visible(False)

    cb_ax = fig.add_axes([0.3, 0.4, 0.45, 0.01])
    sm = plt.cm.ScalarMappable(cmap=pyplot_cmap)
    cbar = fig.colorbar(sm, cax=cb_ax, orientation="horizontal")
    cbar.set_label(label="Metabolic Accesibility", fontsize=14, fontweight="bold")
    cbar.ax.set_xticks(ticks=[0, 0.5, 1.0], labels=["low", "mid", "high"])
    cbar.ax.tick_params(labelsize=12)

    return fig


def display_molecules(data_frame: pd.DataFrame,legends_col: str = "metabolite_id",smiles_col: str = "main_product_smiles",scale_from_column: str = "confidence_score",columns_to_display: list = ["reaction_id"],cmap: str = "RdYlGn"):
    """
    Displays a grid of molecules from a data frame, using different colors to indicate the values of a specified column.

    Parameters
    ----------
    data_frame : pd.DataFrame
        A pandas data frame containing the molecular data.
    legends_col : str, optional
        The name of the column to use as the legend for each molecule. Default is ‘metabolite_id’.
    smiles_col : str, optional
        The name of the column containing the SMILES strings for each molecule. Default is ‘main_product_smiles’.
    scale_from_column : str, optional
        The name of the column to use for scaling the colors of the molecules. Default is ‘confidence_score’.
    columns_to_display : list, optional
        A list of column names to display as tooltips when hovering over the molecules. Default is [‘reaction_id’].
    cmap : str, optional
        The name of the matplotlib colormap to use for coloring the molecules. Default is ‘RdYlGn’.

    Returns
    -------
    mols2grid.display
        A mols2grid display object that shows the grid of molecules with legends, colors and tooltips. The display object has the following features:
            - Each molecule is drawn in a 2D projection with the atom symbols and bond types shown.
            - The legend of each molecule is displayed below the image, using the value from the legends_col column.
            - The color of each molecule is determined by the value from the scale_from_column column, using the cmap colormap. A color bar is displayed on the top right corner of the grid, showing the range of values.
            - The tooltip of each molecule is displayed when hovering over the image, showing the values from the columns_to_display list.
            - The grid can be filtered, sorted and searched by using the widgets on the top left corner of the grid.
    """

    def _cmap_map(function, pyplot_cmap):
        """Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
        This routine will break any discontinuous points in a colormap.

        https://scipy-cookbook.readthedocs.io/items/Matplotlib_ColormapTransformations.html
        """
        cdict = pyplot_cmap._segmentdata
        step_dict = {}
        # Firt get the list of points where the segments start or end
        for key in ("red", "green", "blue"):
            step_dict[key] = list(map(lambda x: x[0], cdict[key]))
        step_list = sum(step_dict.values(), [])
        step_list = np.array(list(set(step_list)))
        # Then compute the LUT, and apply the function to the LUT
        reduced_cmap = lambda step: np.array(pyplot_cmap(step)[0:3])
        old_LUT = np.array(list(map(reduced_cmap, step_list)))
        new_LUT = np.array(list(map(function, old_LUT)))
        # Now try to make a minimal segment definition of the new LUT
        cdict = {}
        for i, key in enumerate(["red", "green", "blue"]):
            this_cdict = {}
            for j, step in enumerate(step_list):
                if step in step_dict[key]:
                    this_cdict[step] = new_LUT[j, i]
                elif new_LUT[j, i] != old_LUT[j, i]:
                    this_cdict[step] = new_LUT[j, i]
            colorvector = list(map(lambda x: x + (x[1],), this_cdict.items()))
            colorvector.sort()
            cdict[key] = colorvector

        return colors.LinearSegmentedColormap("colormap", cdict, 1024)

    color_map = _cmap_map(lambda x: x / 2 + 0.5, plt.get_cmap(cmap))

    if scale_from_column == "confidence_score":
        divnorm = colors.TwoSlopeNorm(vmin=0, vcenter=1.5, vmax=3)
    else:
        divnorm = colors.TwoSlopeNorm(
            vmin=data_frame[scale_from_column].min(),
            vcenter=data_frame[scale_from_column].mean(),
            vmax=data_frame[scale_from_column].max(),
        )

    view = mols2grid.display(
        data_frame,
        smiles_col=smiles_col,
        n_rows=5,
        n_cols=7,
        subset=[legends_col, "img", scale_from_column],
        tooltip=columns_to_display,
        style={
            "__all__": lambda x: f"background-color:{colors.rgb2hex(color_map(divnorm(x[scale_from_column])))};font-weight:bold"
        },
        fixedBondLength=25,
        clearBackground=False,
        gap=1,
    )
    return view


def plot_relationships(data_frame: pd.DataFrame, nodes: list = ["reaction_id", "metabolite_id"]):
    """
    Creates a Sankey diagram to visualize the evidences of metabolite annotations in a data frame.

    Parameters
    ----------
    data_frame : pd.DataFrame
        A pandas data frame that contains the metabolite annotations and their evidences.
    nodes : list, optional
        A list of column names that represent the nodes of the Sankey diagram. The default value is [‘reaction_id’, ‘metabolite_id’].

    Returns
    -------
    plotly.Figure
        A plotly figure object that contains the Sankey diagram. The diagram has the following features:
            - The nodes are arranged horizontally from left to right, corresponding to the order of the columns in the nodes list.
            - The links are drawn as curved lines connecting the nodes, representing the flow of evidences from one node to another.
            - The width of each link is proportional to the number of evidences for that pair of nodes.
            - The color of each link is determined by the color of the source node, using a distinct color for each node.
            - The label of each node is displayed on top of the node, using the value from the corresponding column in the data frame.
            - The tooltip of each link shows the source and target node names and the number of evidences for that link.
            - The title of the diagram shows the names of the columns used for the nodes.
    """
    table = copy.deepcopy(data_frame)
    table = data_frame[nodes]
    for col in table.columns:
        table[col].fillna(f"no_{col}", inplace=True)
        table[col] = table[col].astype(str)
        table[col] = table[col].str.split(";")
        table = table.explode(col, ignore_index=True)

    def get_pairs(input_list):
        result = []
        for i in range(len(input_list) - 1):
            result.append([input_list[i], input_list[i + 1]])
        return result

    def hex_to_rgba(hex_string, alpha=1.0):
        """Converts a hex color value to an RGBA tuple."""
        hex_value = hex_string.lstrip("#")
        rgba = tuple(int(hex_value[i : i + 2], 16) for i in (0, 2, 4)) + (alpha,)
        return rgba

    cols = table.columns
    column_pairs = get_pairs(cols)
    unique_elements = []
    unique_elements += [i for c in cols for i in table[c].unique()]
    indexes = {value: index for index, value in enumerate(unique_elements)}

    sankey_table = pd.DataFrame()
    for i, col_pair in enumerate(column_pairs):
        grouped = table[col_pair].value_counts().reset_index(name="Values")
        grouped.columns = ["Source_name", "Target_name", "Values"]
        sankey_table = pd.concat([sankey_table, grouped], ignore_index=True)

    sankey_table["Source"] = sankey_table.Source_name.map(lambda x: indexes.get(x))
    sankey_table["Target"] = sankey_table.Target_name.map(lambda x: indexes.get(x))

    colors = [
        distinctipy.get_hex(c)
        for c in distinctipy.get_colors(
            n_colors=len(indexes.keys()),
            pastel_factor=0.9,
        )
    ]

    colorsNode = {k: color for k, color in zip(indexes.keys(), colors)}

    colorsNode.update(
        {
            k: colorsNode[table[table[cols[-1]] == k][cols[-2]].iloc[0]]
            for k in table[cols[-1]]
        }
    )
    colorsNode.update(
        {k: "rgba(179, 179, 179,0.3)" for k in table[table.columns[-1]].unique()}
    )

    colorsLinks = []
    for i in sankey_table.index:
        colorsLinks.append(colorsNode[sankey_table.Source_name[i]])

    rgba = [f"rgba{hex_to_rgba(c,0.2)}" for c in colorsLinks]

    data_trace = dict(
        type="sankey",
        domain=dict(x=[0, 1], y=[0, 1]),
        orientation="h",
        valueformat=".0f",
        node=dict(
            pad=20,
            thickness=10,
            line=dict(
                # color = 'red',
                width=1
            ),
            label=list(indexes.keys()),
            color=list(colorsNode.values()),
            # x=x,
            # y=y,
        ),
        link=dict(
            source=sankey_table.Source.to_list(),
            target=sankey_table.Target.to_list(),
            value=sankey_table.Values.to_list(),
            color=rgba,
        ),
    )

    layout = dict(
        title="Evidences",
        title_x=0.5,
        height=600,
        width=1000,
        font=dict(size=12),
    )

    fig = go.Figure(data=data_trace, layout=layout)

    return fig
