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
warnings.filterwarnings('ignore')


class Visualizer:
    """
    Class to plot a molecule with highlighted atoms.
    """
    def __init__(self, data_frame: pd.DataFrame, query: Chem.Mol, name: str = 'Query'):
        self.data_frame=data_frame
        self.query=query
        self.name=name
        
    def plot_confidence_scores(self,cmap: str = 'RdYlGn'):
        fig = px.scatter_3d(self.data_frame, x='similarity_substrates', y='similarity_products', z='reacting_atoms_efficiency',color='confidence_score',range_color=[0,3],
                    color_continuous_scale=cmap,width=1000,height=800,opacity=0.0,hover_name='metabolite_id')

        fig.update_traces(marker_size = 5, showlegend=False,)

        fig.update_layout(paper_bgcolor='white',plot_bgcolor='white',scene = dict(
            xaxis = dict(title = 'Similarity substrates', showgrid = True, showline = True,linecolor='black',backgroundcolor="rgb(200, 200, 230)",gridcolor="white",showbackground=True,zerolinecolor="black", range=(0, 1),dtick=0.25,nticks=4),
            yaxis = dict(title = 'Similarity products', showgrid = True, showline = True,linecolor='black',backgroundcolor="rgb(230, 200,230)",gridcolor="white",showbackground=True,zerolinecolor="black",range=(0, 1),dtick=0.25,nticks=4),
            zaxis = dict(title = 'Atoms efficiency', showgrid = True, showline = True,linecolor='black',backgroundcolor="rgb(230, 230,200)",gridcolor="white",showbackground=True,zerolinecolor="black",range=(0, 1),dtick=0.25,nticks=4)),
                          coloraxis_colorbar_title_text = 'Confidence score', coloraxis_colorbar_len=0.5,coloraxis_colorbar_dtick=0.5
                         )
        return fig

    def plot_evidences(self,nodes:list=['ec','metabolite_id']):
        table=copy.deepcopy(self.data_frame)
        table=self.data_frame[nodes]
        for col in table.columns:
            table[col].fillna(f"no_{col}",inplace=True)
            table[col]=table[col].astype(str)
            table[col]=table[col].str.split(";")
            table=table.explode(col,ignore_index=True)
        
        
        def get_pairs(input_list):
            result = []
            for i in range(len(input_list)-1):
                result.append([input_list[i], input_list[i+1]])
            return result
        
        def hex_to_rgba(hex_string, alpha=1.0):
            """Converts a hex color value to an RGBA tuple."""
            hex_value = hex_string.lstrip('#')
            rgba = tuple(int(hex_value[i:i+2], 16) for i in (0, 2, 4)) + (alpha,)
            return rgba
        
        cols=table.columns
        column_pairs=get_pairs(cols)
        unique_elements=[]
        unique_elements+=[i for c in cols for i in table[c].unique()]
        indexes={value: index for index,value in enumerate(unique_elements)}
        
        self.sankey_table=pd.DataFrame()
        for i,col_pair in enumerate(column_pairs):
            grouped=table[col_pair].value_counts().reset_index(name='Values')
            grouped.columns=['Source_name','Target_name','Values']
            self.sankey_table=self.sankey_table.append(grouped,ignore_index=True)
        
        self.sankey_table['Source']=self.sankey_table.Source_name.map(lambda x: indexes.get(x))
        self.sankey_table['Target']=self.sankey_table.Target_name.map(lambda x: indexes.get(x))
        
        colors=[distinctipy.get_hex(c) for c in distinctipy.get_colors(n_colors=len(indexes.keys()),pastel_factor=0.9,)]
        
        colorsNode={k:color for k,color in zip(indexes.keys(),colors)}
        
        colorsNode.update({k:colorsNode[table[table[cols[-1]]==k][cols[-2]].iloc[0]] for k in table[cols[-1]]})
        colorsNode.update({k:'rgba(179, 179, 179,0.3)' for k in table[table.columns[-1]].unique()})
        
        colorsLinks=[]
        for i in self.sankey_table.index:
            colorsLinks.append(colorsNode[self.sankey_table.Source_name[i]])

        rgba=[f"rgba{hex_to_rgba(c,0.2)}" for c in colorsLinks]

        data_trace = dict(type='sankey',domain = dict(x = [0,1],y =  [0,1]),
                              orientation = "h", valueformat = ".0f",
                              node = dict(pad = 20,
                                          thickness = 10,
                                          line = dict(
                                              #color = 'red',
                                              width = 1),
                                          label =  list(indexes.keys()),
                                          color = list(colorsNode.values()),
                                          #x=x,
                                          #y=y,
                                         ),
                              link = dict(
                                  source = self.sankey_table.Source.to_list(),
                                  target = self.sankey_table.Target.to_list(),
                                  value = self.sankey_table.Values.to_list(),
                                  color = rgba,
                              ))

        layout = dict(
            title = "Evidences",title_x=0.5,
            height = 800,width=1400,
            font = dict(size =12),)

        fig = go.Figure(data=data_trace, layout=layout)
        
        return fig
    
    def plot_metabolic_accesibility(self,cmap: str = 'RdYlGn_r'):
        reacting_atoms=[a for l in self.data_frame.reacting_atoms_in_query for a in l]
        pyplot_cmap = plt.get_cmap(cmap)
        atoms_counts = {i: reacting_atoms.count(i) for i in reacting_atoms}
        
        def _scale_value(val, values):
            """
            Scales the value between 0 and 1.
            """
            return (val - min(values)) / (max(values) - min(values))
        
        self.color_atom_map = {k: [pyplot_cmap(_scale_value(v, atoms_counts.values()), alpha=0.1)] for k, v in atoms_counts.items()}
        d2d = rdMolDraw2D.MolDraw2DCairo(1600,400)
        d2d.ClearDrawing()
        dos = d2d.drawOptions()
        dos.bondLineWidth = 8
        dos.centreMoleculesBeforeDrawing = True
        dos.highlightRadius = 0.5
        dos.useBWAtomPalette()
        d2d.DrawMoleculeWithHighlights(
            self.query, 
            legend='', 
            highlight_atom_map=self.color_atom_map, 
            highlight_bond_map={}, 
            highlight_radii={}, 
            highlight_linewidth_multipliers={}
        )

        d2d.FinishDrawing()
        img = Image.open(io.BytesIO(d2d.GetDrawingText()))

        fig, axes = plt.subplots(figsize=(8,15))
        ax = plt.imshow(img)
        ax.axes.set_title(self.name, fontsize=18, fontweight='bold')
        ax.axes.set_xticklabels([])
        ax.axes.set_yticklabels([])
        ax.axes.tick_params(axis='both', which='both', length=0)
        for spine in ['bottom','top','right','left']:
            ax.axes.spines[spine].set_visible(False)

        cb_ax = fig.add_axes([0.3, 0.4, 0.45, 0.01])
        sm = plt.cm.ScalarMappable(cmap=pyplot_cmap)
        cbar = fig.colorbar(sm, cax=cb_ax, orientation='horizontal')
        cbar.set_label(label='Metabolic Accesibility', fontsize=14, fontweight='bold')
        cbar.ax.set_xticks(ticks=[0,0.5,1.0], labels=['low','mid','high'])
        cbar.ax.tick_params(labelsize=12)

        return fig
        
    def display_metabolites(self, legends_col:str = 'metabolite_id', smiles_col:str = 'main_product_smiles',scale_from_column:str = 'confidence_score',columns_to_display:list = ['mnx_id'],cmap:str='RdYlGn'):        
        
        def _cmap_map(function, pyplot_cmap): 
            """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
            This routine will break any discontinuous points in a colormap.

            https://scipy-cookbook.readthedocs.io/items/Matplotlib_ColormapTransformations.html
            """
            cdict = pyplot_cmap._segmentdata
            step_dict = {}
            # Firt get the list of points where the segments start or end
            for key in ('red', 'green', 'blue'):
                step_dict[key] = list(map(lambda x: x[0], cdict[key]))
            step_list = sum(step_dict.values(), [])
            step_list = np.array(list(set(step_list)))
            # Then compute the LUT, and apply the function to the LUT
            reduced_cmap = lambda step : np.array(pyplot_cmap(step)[0:3])
            old_LUT = np.array(list(map(reduced_cmap, step_list)))
            new_LUT = np.array(list(map(function, old_LUT)))
            # Now try to make a minimal segment definition of the new LUT
            cdict = {}
            for i, key in enumerate(['red','green','blue']):
                this_cdict = {}
                for j, step in enumerate(step_list):
                    if step in step_dict[key]:
                        this_cdict[step] = new_LUT[j, i]
                    elif new_LUT[j,i] != old_LUT[j, i]:
                        this_cdict[step] = new_LUT[j, i]
                colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
                colorvector.sort()
                cdict[key] = colorvector

            return colors.LinearSegmentedColormap('colormap',cdict,1024)
        
        color_map = _cmap_map(lambda x: x/2 + 0.5,plt.get_cmap(cmap))
        
        if scale_from_column=='confidence_score':
            divnorm = colors.TwoSlopeNorm(vmin=0, vcenter=1.5, vmax=3)
        else:
            divnorm = colors.TwoSlopeNorm(vmin=self.data_frame[scale_from_column].min(), vcenter=self.data_frame[scale_from_column].mean(), vmax=self.data_frame[scale_from_column].max())
        
        view = mols2grid.display(self.data_frame, smiles_col=smiles_col, n_rows=5, n_cols=7,
                          subset=[legends_col, "img", scale_from_column],
                          tooltip=columns_to_display,
                          style={
                                "__all__": lambda x: f"background-color:{colors.rgb2hex(color_map(divnorm(x[scale_from_column])))};font-weight:bold"
                          },
                          fixedBondLength=25, clearBackground=False, gap=1)
        
        return view
    