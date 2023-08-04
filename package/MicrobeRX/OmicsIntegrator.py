"""
Functionalities to analyze the enzymes and organisms involved in metabolism of Xenobiotics
====================================================
"""

import plotly.express as px
from plotly.subplots import make_subplots
from distinctipy import distinctipy

import pandas as pd

from package.MicrobeRX.DataFiles import load_microbes_data, load_microbes_reactions


__all__=['plot_species_sunburst']
''' 
################### TOOLKIT ###################
'''

def check_if_microbes_databases_are_loaded():
    
    global MICROBES_DATA
    MICROBES_DATA=None
    if MICROBES_DATA is not None:
        pass
    else:       
        MICROBES_DATA=load_microbes_data()

    global MICROBES_REACTIONS
    MICROBES_REACTIONS=None
    if MICROBES_REACTIONS is not None:
        pass
    else:
        MICROBES_REACTIONS=load_microbes_reactions()



def plot_species_sunburst(sources:list,path:str='short'):
    
    check_if_microbes_databases_are_loaded()
    
    selection=MICROBES_REACTIONS[sources].dropna(subset=sources,axis=0,how='all')
    
    subplots = len(selection.columns)
    cols = 3
    rows = subplots // cols 
    if subplots-(rows*cols) >0:
        rows += 1

    matrix = [dict(row=r+1,column=c+1) for r in range(rows) for c in range(cols)]

    specs=[[{"type": "domain"} for c in range(cols)] for r in range(rows)]

    F = make_subplots(rows=rows, cols=cols,horizontal_spacing=0.02,vertical_spacing=0.1,specs=specs,subplot_titles=selection.columns.to_list())
    
    phyla=MICROBES_DATA.Phylum.unique()
    colors=[distinctipy.get_hex(c) for c in distinctipy.get_colors(n_colors=len(phyla),pastel_factor=0.8)]
    color_map=dict(zip(phyla,colors))
    
    if path == 'full':
        path=['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    if path == 'short':
        path=['Kingdom', 'Phylum', 'Order', 'Genus', 'Species']
    
    for index,col in enumerate(selection.columns):
        slide=MICROBES_DATA[MICROBES_DATA.microbe_name.isin(selection[[col]].dropna().index)]
        
        fig1=px.sunburst(slide,path=path,color='Phylum',color_discrete_map=color_map,branchvalues='total')

        for d in fig1.data:
            F.add_trace(d,row=matrix[index]['row'],col=matrix[index]['column'])

    F.update_layout(grid= dict(rows=rows,columns=cols,),height=400*rows, width=400*cols, title_text=f"Species Sunburst")

    return F