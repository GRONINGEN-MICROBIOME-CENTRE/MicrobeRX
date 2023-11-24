"""
This is a module that provides functions to analyze organims, and enzyme involved in the production of metabolites predicted by MicrobeRX.

The module contains the following functions:

- plot_species_sunburst: This function creates a sunburst plot of the microbial species in the sources list. 

- fetch_batch_sequences: This function fetches a list of sequences from the NCBI Entrez database. It uses the Biopython library to access the Entrez API and parse the FASTA format. It also uses a helper function _fetch_sequence to fetch and return a single sequence.

- get_interpro: This function retrieves the InterProScan results for a given sequence from the EBI InterProScan 5 web service. 

- plot_interpro_results: This functio creates a bar plot of the InterProScan results for a given sequence.

- run_multi_sequence_aligment: This function performs a multiple sequence alignment (MSA) and a phylogenetic tree construction for a given set of sequences using the ClustalW2 program.

- plot_similarity_matrix: This function creates a heatmap of the pairwise similarity scores for a given set of sequences using the Dash Bio library. It also accepts optional parameters to choose the color map and the homology percentage for the heatmap.

- plot_aligment_chart: This function creates a chart of the multiple sequence alignment (MSA) for a given set of sequences using the Dash Bio library. It also accepts optional parameters to choose the color scale and the conservation method for the chart.

"""

__all__ = ["plot_species_sunburst","fetch_batch_sequences","get_interpro","plot_interpro_results","run_multi_sequence_aligment","plot_similarity_matrix","plot_aligment_chart"]

import requests, time, re, copy, io

from importlib_resources import files

import plotly.express as px
from plotly.subplots import make_subplots
import dash_bio
from dash import html
from jupyter_dash import JupyterDash

from distinctipy import distinctipy

import pandas as pd

from Bio import Entrez
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO, SeqRecord

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

    return selection, F



def fetch_batch_sequences(entries:list=None, sequence_ids:list=None, email:str=None, database:str="protein"):
    """ 
    The function fetch_batch_sequences fetches a list of sequences from the NCBI Entrez database. It uses the Biopython library to access the Entrez API and parse the FASTA format. It also uses a helper function _fetch_sequence to fetch and return a single sequence.

    Parameters
    ----------
        entries: A list of strings that represent the accession numbers of the sequences to be fetched. For example, ['WP_015582217.1', 'WP_001277567.1'].
        sequence_ids: A list of strings that represent the custom ids to be assigned to the fetched sequences. For example, ['seq1', 'seq2']. The length of this list should match the length of the entries list.
        email: A string that specifies the email address of the user. This is required by the Entrez API to identify the user and avoid abusing the system. For example, 'user@example.com'.
        database: A string that specifies the name of the Entrez database to fetch the sequences from. The default value is 'protein'. For example, 'nucleotide'.

    Returns
    -------
    sequences: A list of Bio.SeqRecord objects that contain the fetched sequences. Each sequence has the id attribute set to the corresponding value in the sequence_ids list. If an error occurs while fetching a sequence, it is skipped and not added to the list.
    """
    
    def _fetch_sequence(entry:str=None, sequence_id:str=None, email:str=None, database:str="protein"):
        Entrez.email = email
        handle = Entrez.efetch(db=database, id=entry, rettype="fasta", retmode="text")
        sequence=SeqIO.read(io.StringIO(handle.read()),format='fasta')
        sequence.id=sequence_id
        handle.close()

        return sequence
    
    sequences=[]
    for index,entry in enumerate(entries):
        try:
            sequence=_fetch_sequence(entry=entries[index],sequence_id=sequence_ids[index],email=email)
            sequences.append(sequence)
        except Exception:
            pass
    
    return sequences

def get_interpro(sequence_id:str=None, sequence:str=None,email:str=None, sequence_type:str='protein',go_terms:bool=True, pathways:bool=True):
    """ 
    The function get_interpro retrieves the InterProScan results for a given sequence from the EBI InterProScan 5 web service. It uses the requests library to access the REST API and the pandas library to parse the tab-separated values (TSV) format. It also accepts optional parameters to include GO terms and pathway information in the output.

    Parameters
    ----------
        sequence_id: A string that represents the id of the sequence to be scanned. For example, 'seq1'.
        sequence: A string of the sequence to be scanned. For example, MKKLLIISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCP.
        email: A string that specifies the email address of the user. This is required by the EBI InterProScan 5 web service to identify the user and avoid abusing the system. For example, 'user@example.com'.
        sequence_type: A string that specifies the type of the sequence to be scanned. It can be either 'protein' or 'nucleotide'. The default value is 'protein'.
        go_terms: A boolean that indicates whether to include GO terms in the output. The default value is True.
        pathways: A boolean that indicates whether to include pathway information in the output. The default value is True.

    Returns
    -------
    interpro: A pandas DataFrame object that contains the InterProScan results. It has the following columns: ['accesion', 'token', 'sequence_length', 'analysis', 'signature_accession', 'signature_description', 'start_location', 'stop_location', 'score', 'status', 'date', 'interpro_accession', 'interpro_description', 'go_annotations', 'pathways']. The last two columns are optional depending on the values of the go_terms and pathways parameters. If an error occurs while fetching the results, it returns the error message as a string.
    """

    if sequence_type=="protein":
        sequence_type="p"
    elif sequence_type=="nucleotide":
        sequence_type="n"
    else:
        print("Sequence must be 'protein' or 'nucleotide'")
    
    
    columns=['accesion','token','sequence_length','analysis','signature_accession','signature_description','start_location','stop_location','score','status','date','interpro_accession','interpro_description']
    
    
    if go_terms == True:
        go_terms="true"
        columns.append('go_annotations')
    elif go_terms == False:
        go_terms="false"
    else:
        print("Go Terms parameter must be bool type True or False")
    
    if pathways == True:
        pathways="true"
        columns.append('pathways')
    elif pathways == False:
        pathways="false"
    else:
        print("Pathways parameter must be bool type True or False")
    
    
    params = {
        "email" :email,
        "stype":sequence_type,
        "goterms": go_terms, # Include GO terms in the output
        "pathways": pathways, # Include pathway information in the output
        "sequence": f">{sequence_id}\n{sequence}", # The protein sequence in fasta format
    }
    
    
    base_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/"

    response = requests.post(base_url +"run", data=params)

    if response.status_code == 200:

        job_id = response.text
        print("Job ID:", job_id)

        while True:

            status = requests.get(base_url + "status/" + job_id).text
            if status == "QUEUED":
                print("Job is queued, please wait...")
                time.sleep(10) 

            elif status == "RUNNING":
                print("Job is running, please wait...")
                time.sleep(30) 
            elif status == "FINISHED":
                results = requests.get(base_url + "result/" + job_id + "/tsv").text
                interpro=pd.read_csv(io.StringIO(results),sep='\t',names=columns)
                break 
            else:
                print("Job failed or not found, status:", status)
                break
    else:
        return response.text
    
    return interpro

def plot_interpro_results(interpro_results:pd.DataFrame=None, compact:str=True):
    """ 
    The function plot_interpro_results creates a bar plot of the InterProScan results for a given sequence. It uses the Plotly Express library to create the bar plot. It also accepts an optional parameter to choose between a compact or a detailed view of the results.

    Parameters
    ----------
        interpro_results: A pandas DataFrame object that contains the InterProScan results. It should have the following columns: ['accesion', 'token', 'sequence_length', 'analysis', 'signature_accession', 'signature_description', 'start_location', 'stop_location', 'score', 'status', 'date', 'interpro_accession', 'interpro_description', 'go_annotations', 'pathways']. The last two columns are optional depending on the values of the go_terms and pathways parameters in the get_interpro function.
        compact: A boolean that indicates whether to use a compact or a detailed view of the results. The default value is True. If True, the bar plot shows the InterPro accession and description for each segment of the sequence. If False, the bar plot shows the analysis and signature description for each segment of the sequence.

    Returns
    -------
    fig: A Plotly Figure object that contains the bar plot. It has one subplot for the sequence. The bar plot shows the distribution of the InterProScan results by their start and stop locations on the sequence. The color of each segment is determined by the InterPro accession or the analysis depending on the value of the compact parameter. The text of each segment shows the InterPro description or the signature description depending on the value of the compact parameter.
    """
    
    interpro_results['sequence_delta'] = interpro_results.stop_location - interpro_results.start_location
    interpro_results.sort_values(by='sequence_delta',ascending=False,ignore_index=True,inplace=True)
    
    if compact==True:
        interpro_compact=copy.deepcopy(interpro_results)
        interpro_compact['all_analysis']=interpro_compact.interpro_accession.map(interpro_compact.groupby('interpro_accession')['analysis'].agg(lambda x: '; '.join(i for i in set(x))))
        interpro_compact.drop_duplicates(subset='interpro_accession',inplace=True,ignore_index=True)
        interpro_compact=interpro_compact[interpro_compact.interpro_accession!='-']
        
        fig = px.bar(interpro_compact, 
             base = "start_location",
             x = "sequence_delta",
             y = "all_analysis",
             color = "interpro_accession", 
             orientation = 'h', 
             text='interpro_description',
             color_discrete_sequence=px.colors.qualitative.Pastel,
             opacity=0.4,
            )


        fig.update_yaxes(autorange="reversed",title = "")
        fig.update_xaxes(title = 'Sequence', showline = True, linecolor='darkgray',linewidth=2.5,range=[1, interpro_compact.sequence_length.iloc[0]])
        fig.update_traces(textposition='inside',textfont_size=12)
        fig.update_layout(title=interpro_compact.accesion.iloc[0],
                          bargap=0.1,
                          height=350,
                          width=1000,
                          paper_bgcolor="white",
                          plot_bgcolor="white",
                          legend=dict(title='InterPro')
                         )
        return fig
    
    if compact==False:
        fig = px.bar(interpro_results, 
             base = "start_location",
             x = "sequence_delta",
             y = "analysis",
             color = "analysis", 
             orientation = 'h', 
             text='signature_description',
             color_discrete_sequence=px.colors.qualitative.Pastel,
             opacity=0.4,
            )


        fig.update_yaxes(autorange="reversed",title = "")
        fig.update_xaxes(title = 'Sequence', showline = True, linecolor='darkgray',linewidth=2.5,range=[1, interpro_results.sequence_length.iloc[0]])
        fig.update_traces(textposition='inside',textfont_size=12)
        fig.update_layout(title=interpro_results.accesion.iloc[0],
                          bargap=0.1,
                          height=350,
                          width=1000,
                          paper_bgcolor="white",
                          plot_bgcolor="white",
                          legend=dict(title='Software')
                         )
        return fig
    

def run_multi_sequence_aligment(sequences_file:str=None,input_format:str="fasta",output_aligment_format:str="fasta"):
    """ 
    The function run_multi_sequence_aligment performs a multiple sequence alignment (MSA) and a phylogenetic tree construction for a given set of sequences using the ClustalW2 program. It uses the Biopython library to parse the input and output files and to run the ClustalW2 command line. It also returns a heatmap of the pairwise alignment scores.

    Parameters
    ----------
        sequences_file: A string that represents the name of the file that contains the sequences to be aligned. For example, 'sequences.faa'.
        input_format: A string that specifies the format of the input file. The default value is 'fasta'. For example, 'phylip'.
        output_aligment_format: A string that specifies the format of the output alignment file. The default value is 'fasta'. For example, 'clustal'.
    
    Returns
    -------
    similarity_matrix: A pandas DataFrame object that contains the heatmap of the pairwise alignment scores. It has the sequence ids as the row and column labels. The values are the percentage of identical positions in the pairwise alignment. The diagonal values are 100. For example:

              seq1  seq2  seq3
        seq1  100.0  85.0  75.0
        seq2   85.0 100.0  80.0
        seq3   75.0  80.0 100.0

    Side effects
    ------------
    The function also creates two output files in the same directory as the input file: 
    - A file named 'sequences.fasta' that contains the MSA in the specified output format. For example, 'sequences.fasta'.
    - A file named 'sequences.dnd' that contains the phylogenetic tree in the Newick format. For example:

        (((seq1:0.02941,seq2:0.02941):0.02941,seq3:0.05882):0.00000,);
    """
    
    records={r.id:r for r in SeqIO.parse(sequences_file, input_format)}
    
    CLUSTALW2=files("microberx.bin").joinpath("clustalw2")
    
    clustalw_cline = ClustalwCommandline(CLUSTALW2, infile=sequences_file,output=output_aligment_format,score='PERCENT',outorder='ALIGNED')
    
    raw_results=clustalw_cline()

    aligment_results=raw_results[0].split("\n")
    
    similarity_matrix=pd.DataFrame()
    
    print(f'MSA and Phylogenetic tree have been saven in {sequences_file} directory')
    
    for line in aligment_results:
        if "Aligned. Score" in line:
            index_match=re.findall(r"(\d+):(\d+)",line)
            indexes=[int(i) for i in index_match[0]]

            score_match=re.findall(r"Score:\s+(\d+)",line)
            score=int(score_match[0])

            similarity_matrix.loc[indexes[0],indexes[1]]=score
            similarity_matrix.loc[indexes[1],indexes[0]]=score
    
    similarity_matrix.fillna(value=100,inplace=True)
    
    similarity_matrix.columns=records.keys()
    similarity_matrix.index=records.keys()
    
    return similarity_matrix

def plot_similarity_matrix(similarity_matrix:pd.DataFrame=None,homology_percentage:float=90,cmap:str='custom'):
    """ 
    The function plot_similarity_matrix creates a heatmap of the pairwise similarity scores for a given set of sequences using the Dash Bio library. It also accepts optional parameters to choose the color map and the homology percentage for the heatmap.

    Parameters
    ----------
        similarity_matrix: A pandas DataFrame object that contains the pairwise similarity scores for the sequences. It should have the sequence ids as the row and column labels. The values should be the percentage of identical positions in the pairwise alignment. For example:

              seq1  seq2  seq3
        seq1  100.0  85.0  75.0
        seq2   85.0 100.0  80.0
        seq3   75.0  80.0 100.0

        homology_percentage: A float that specifies the threshold for the homology color. The default value is 90. It should be between 0 and 100. The dendogram will use a different color for the values that are above or equal to this threshold displaying homology clusters.
        cmap: A string or a list that specifies the color map for the heatmap. The default value is 'custom'. If 'custom', the color map is [[0.0, 'rgb(64, 126, 156)'], [0.5,'rgb(242,241,241)'], [1.0, 'rgb(195,85,58)']]. Otherwise, it can be one of the predefined color maps in the Dash Bio library: 'Blackbody', 'Bluered', 'Blues', 'Earth', 'Electric', 'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow', 'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd'. Alternatively, it can be a custom color map as a list of lists that map a value between 0 and 1 to a color. For example, [[0.0, 'red'], [0.5, 'white'], [1.0, 'blue']].

    Returns
    -------
    fig: A Plotly Figure object that contains the heatmap of the pairwise similarity scores. It has the following features:
    - It shows the sequence ids and the similarity scores for each pair of sequences in the heatmap.
    - It allows the user to zoom in and out, pan, and select a region of the heatmap.
    - It allows the user to change the color map, the homology percentage, and the display options of the heatmap.
    - It shows a color bar that indicates the range of the similarity scores and the homology color.
    """
    
    if cmap=="custom":
        cmap=[[0.0, 'rgb(64, 126, 156)'], [0.5,'rgb(242,241,241)'], [1.0, 'rgb(195,85,58)']]
    
    fig=dash_bio.Clustergram(
        data=similarity_matrix.values,
        column_labels=list(similarity_matrix.columns.values),
        row_labels=list(similarity_matrix.index),
        height=round(len(similarity_matrix.index)*20),
        width=900,
        color_map=cmap,
        line_width=1.5,
        hidden_labels='col',
        display_ratio=0.2,
        center_values=False,
        color_threshold={'row':homology_percentage,'col': homology_percentage})

    heatmap_trace=fig.data[-1]
    heatmap_trace.update(zmin=0, zmax=100, colorbar= {'title':"% Identity","x":-0.2,"y":0.95,"orientation":"v","len":0.3,"thickness":25,"dtick":25, "tick0":0, "nticks":4, "tickfont":{"size":10}})

    return fig

def plot_aligment_chart(msa_file:str=None, cmap:str="custom",color_scale:str='mae'):
    
    """ 
    The function plot_aligment_chart creates a chart of the multiple sequence alignment (MSA) for a given set of sequences using the Dash Bio library. It also accepts optional parameters to choose the color scale and the conservation method for the chart.

    Parameters
    ----------
        msa_file: A string that represents the name of the file that contains the MSA in the FASTA format. For example, 'sequences.fasta'.
        cmap: A string or a list that specifies the color map for the conservation scores. The default value is 'custom'. If 'custom', the color map is [[0.0, 'rgb(64, 126, 156)'], [0.5,'rgb(242,241,241)'], [1.0, 'rgb(195,85,58)']]. Otherwise, it can be one of the predefined color maps in the Plotly library: 'viridis', 'RdBu', etc...
        color_scale: A string that specifies the color scale for the alignment symbols. The default value is 'mae'. It can be one of the predefined color scales in the Dash Bio library: 'buried', 'cinema', 'clustal', 'clustal2', 'helix', 'hydrophobicity', 'lesk', 'mae', 'nucleotide', 'purine', 'strand', 'taylor', 'turn', or 'zappo'. Alternatively, it can be a custom color map as a dictionary that maps each nucleotide or amino acid to a color. For example, {'A': 'red', 'C': 'blue', 'G': 'green', 'T': 'yellow'}.

    Returns
    -------
    None

    Side effects
    ------------
    The function also creates and runs a Dash app that displays the chart of the MSA. The chart has the following features:
    - It shows the sequence ids, the alignment symbols, and the conservation scores for each position in the alignment.
    - It allows the user to zoom in and out, pan, and select a region of the alignment.
    """
    
    if cmap=="custom":
        cmap=[[0.0, 'rgb(64, 126, 156)'], [0.5,'rgb(242,241,241)'], [1.0, 'rgb(195,85,58)']]
    
    with open(msa_file,'r') as f:
        data=f.read()
    
    chart=dash_bio.AlignmentChart(id='MicrobeRX',data=data,textsize=8,correctgap=True,
                                  colorscale=color_scale,conservationcolorscale=cmap,
                                  conservationmethod='conservation',numtiles=50,
                                  groupbars=False,showid=False,showgap=True,
                                  width=1000,height=round(45*30))

    app = JupyterDash(__name__)

    app.layout = html.Div(chart)


    app.run()