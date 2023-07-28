import requests, copy, random
import pandas as pd
import pubchempy as pcp
from pyopenms import *

import plotly.express as px
import plotly.graph_objects as go

from distinctipy import distinctipy

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from tqdm.notebook import tqdm


''' 
################### TOOLKIT ###################
'''

def compute_molecular_descriptors(data_frame: pd.DataFrame, smiles_col:str) -> pd.DataFrame:
    ''' 
    This function computes some molecular descriptors for a given data frame of SMILES strings.
    Parameters:
    data_frame: a pandas data frame that contains SMILES strings of molecules
    smiles_col: the name of the column that contains the SMILES strings
    Returns:

    data_frame: the same data frame as input, but with additional columns for each molecular descriptor computed. The descriptors are:
        MolWt: the molecular weight of the molecule
        LogP: the octanol-water partition coefficient of the molecule
        NumHAcceptors: the number of hydrogen bond acceptors in the molecule
        NumHDonors: the number of hydrogen bond donors in the molecule
        NumRotatableBonds: the number of rotatable bonds in the molecule
        TPSA: the topological polar surface area of the molecule
        MolFormula: the molecular formula of the molecule 
    '''
    data_frame=copy.deepcopy(data_frame)
    for index in (data_frame.index):
        mol=Chem.MolFromSmiles(data_frame[smiles_col][index])
        if mol:
            data_frame.loc[index,'MolWt']= round(Descriptors.MolWt(mol),3) 
            data_frame.loc[index,'LogP']= round(Descriptors.MolLogP(mol),3) 
            data_frame.loc[index,'NumHAcceptors']= Descriptors.NumHAcceptors(mol) 
            data_frame.loc[index,'NumHDonors']= Descriptors.NumHDonors(mol) 
            data_frame.loc[index,'NumRotatableBonds']= Descriptors.NumRotatableBonds(mol)
            data_frame.loc[index,'TPSA']= Descriptors.TPSA(mol) 
            data_frame.loc[index,'MolFormula']= AllChem.CalcMolFormula(mol)
    
    return data_frame

def plot_molecular_descriptors(data_frame:pd.DataFrame, names_col:str) -> go.Figure:
    '''
    This function plots the molecular descriptors of a given data frame using polar coordinates.

    Parameters:
    - data_frame: a pandas data frame that contains the molecular descriptors as columns and the compound names as rows.
    - names_col: a string that specifies the name of the column that contains the compound names.

    Returns:
    - Figure: a plotly figure object that shows the polar plot of the molecular descriptors.

    The function normalizes the molecular descriptors to fit in the range [0, 1] and then plots them as radial lines for each compound. The function also plots the upper and lower limits of the Lipinski's rule of five as shaded regions in orange and yellow, respectively. The function uses distinct colors for each compound and displays a legend on the right side of the plot.
    '''
        
    data=pd.DataFrame(index=list(data_frame[names_col])) 
    data['MolWt']=[i/500 for i in data_frame['MolWt']]
    data['LogP']=[i/5 for i in data_frame['LogP']]
    data['HBA']=[i/10 for i in data_frame['NumHAcceptors']]
    data['HBD']=[i/5 for i in data_frame['NumHDonors']]
    data['RotB']=[i/10 for i in data_frame['NumRotatableBonds']]
    data['TPSA']=[i/140 for i in data_frame['TPSA']]

    Ro5_up=[1,1,1,1,1,1]
    Ro5_low=[0.5,0.1,0.1,0.25,0.1,0.5]  
    #data=data.reindex(natsorted(data.index))
    categories=list(data.columns)  

    Figure = go.Figure()

    fig1 = px.line_polar(r=Ro5_up, theta=categories, line_close=True,color_discrete_sequence=['red'],)
    fig1.update_traces(fill='toself',fillcolor='orange',opacity=0.4,name='Upper limit')
    fig2 = px.line_polar(r=Ro5_low, theta=categories, line_close=True,color_discrete_sequence=['red'])
    fig2.update_traces(fill='toself',fillcolor='yellow',opacity=0.4,name='Lower limit')

    fig1.data[0].showlegend=True
    fig2.data[0].showlegend=True

    Figure.add_trace(fig1.data[0])
    Figure.add_trace(fig2.data[0])


    color = color = [distinctipy.get_hex(c) for c in distinctipy.get_colors(n_colors=len(data.index),pastel_factor=0.7)]

    for i,name in enumerate(data.index):
        f=px.line_polar(r=list(data.loc[name]), theta=categories, line_close=True,color_discrete_sequence=[color[i]])
        f.update_traces(name=name)
        f.data[0].showlegend=True
        Figure.add_trace(f.data[0])

    Figure.update_layout( width=1200,height=700,
      polar=dict(
        radialaxis=dict(
          visible=True)),showlegend=True,legend=dict(x=1.2, y=0.95),legend_orientation="h")

    return Figure

def compute_isotopic_mass(data_frame:pd.DataFrame, molformula_col:str) -> pd.DataFrame:
    '''
    This function computes the isotopic mass distribution of a given data frame using the pyOpenMS library.

    Parameters:
    - data_frame: a pandas data frame that contains the molecular formulas as a column.
    - molformula_col: a string that specifies the name of the column that contains the molecular formulas.

    Returns:
    - data_frame: a pandas data frame that has two additional columns: 'probability_sum' and 'mass_distribution'.

    The function iterates over the rows of the data frame and uses the EmpiricalFormula class from pyOpenMS to create an object for each molecular formula. Then, it uses the CoarseIsotopePatternGenerator class to generate the isotopic mass distribution with a resolution of 4. It calculates the sum of the probabilities of all isotopes and stores it in the 'probability_sum' column. It also formats the mass and probability of each isotope as a string and stores it in the 'mass_distribution' column, separated by semicolons. If an exception occurs, the function skips the row and continues with the next one.
    
    Example:
    >>> import pandas as pd
    >>> from pyopenms import EmpiricalFormula, CoarseIsotopePatternGenerator
    >>> df = pd.DataFrame({'formula': ['C6H12O6', 'C2H4O2', 'C3H8O3']})
    >>> df = compute_isotopic_mass(df, 'formula')
    >>> print(df)
        formula  probability_sum                                   mass_distribution
    0  C6H12O6            1.0000  180.0634:100.0;181.0668:10.72;182.0701:1.176;183...
    1   C2H4O2            1.0000  60.0211:100.0;61.0245:11.08;62.0279:1.216;63.031...
    2   C3H8O3            0.9999  92.0473:100.0;93.0507:10.55;94.0541:1.159;95.057...
    '''
    data_frame=copy.deepcopy(data_frame)      
    for index in (data_frame.index):
        try:
            molF = EmpiricalFormula(data_frame[molformula_col][index])
            isotopes=molF.getIsotopeDistribution(CoarseIsotopePatternGenerator(4))
            prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
            distribution=';'.join([f'{round(iso.getMZ(),4)}:{round(iso.getIntensity()*100,4)}' for iso in isotopes.getContainer()])
            sumR=round(prob_sum,4)  
            data_frame.loc[index,'probability_sum']=sumR
            data_frame.loc[index,'mass_distribution']=distribution
        except Exception:
            pass
    
    return data_frame

def plot_isotopic_masses(data_frame:pd.DataFrame, names_col:str,mass_distribution_col:str) -> go.Figure:
    '''
    This function plots the isotopic mass distribution of a given data frame using plotly.

    Parameters:
    - data_frame: a pandas data frame that contains the isotopic mass distribution as a column of strings, where each string has the format 'mass:probability;mass:probability;...'
    - names_col: a string that specifies the name of the column that contains the compound names.
    - mass_distribution_col: a string that specifies the name of the column that contains the isotopic mass distribution.

    Returns:
    - Figure: a plotly figure object that shows the bar plot of the isotopic mass distribution for each compound.

    The function creates a copy of the data frame and drops any rows that have missing values in the mass distribution column. Then, it converts the mass distribution column from strings to dictionaries, where the keys are the masses and the values are the probabilities. It also resets the index of the data frame. Then, it iterates over the rows of the data frame and creates a bar plot for each compound using plotly express. It uses distinct colors for each compound and displays a legend on the right side of the plot.
    '''
    
    masses=data_frame[[names_col,mass_distribution_col]]
    masses.dropna(subset=[mass_distribution_col],inplace=True)
    masses.mass_distribution=masses[mass_distribution_col].apply(lambda x: {float(value.split(':')[0]):float(value.split(':')[1]) for value in x.split(';')})
    masses.reset_index(drop=True,inplace=True)

    Figure = go.Figure()

    color = [distinctipy.get_hex(c) for c in distinctipy.get_colors(n_colors=len(masses.index),pastel_factor=0.7)]

    for index in masses.index:
        f=px.bar(x=tuple(masses[mass_distribution_col][index].keys()), y=tuple(masses[mass_distribution_col][index].values()),color_discrete_sequence=[color[index]])
        f.update_traces(name=masses[names_col][index])
        f.data[0].showlegend=True
        Figure.add_trace(f.data[0])

    Figure.update_layout(width=1200,height=700,showlegend=True,legend=dict(x=1.2, y=0.95),legend_orientation="h",plot_bgcolor='rgba(0, 0, 0, 0)',font=dict(size=12),xaxis_title="Atomic Mass (u)",yaxis_title="Relative abundance (%)",title="Isotopic distribution")

    return Figure

def search_pubchem(data_frame,entry_col:str,entry_type:str='smiles') -> pd.DataFrame:
    '''
    This function searches the PubChem database for compounds that match a given data frame of identifiers.

    Parameters:
    - data_frame: a pandas data frame that contains the identifiers of the compounds to search for.
    - entry_col: a string that specifies the name of the column that contains the identifiers.
    - entry_type: a string that specifies the type of the identifiers, such as 'smiles', 'inchi', 'cid', etc. The default is 'smiles'.

    Returns:
    - data_frame: the same data frame as input, but with additional columns for each PubChem property retrieved. The properties are:
        PubChem_CID: the PubChem compound identifier, separated by semicolons if there are multiple matches.
        PubChem_SID: the PubChem substance identifier, separated by semicolons if there are multiple matches. Only the first three SIDs are shown.
        PubChem_Synonyms: the synonyms of the compound, separated by semicolons if there are multiple matches.

    The function iterates over the rows of the data frame and uses the pubchempy library to query the PubChem database for compounds that match the identifier in the specified column and namespace. It extracts the CIDs, SIDs and synonyms of the matching compounds and stores them in the corresponding columns of the data frame. If an exception occurs, the function skips the row and continues with the next one.
    '''
    for index in tqdm(data_frame.index):
        try:
            matches=pcp.get_compounds(data_frame[entry_col][index],namespace=entry_type)
            cids=[mat.cid for mat in matches]
            sids=[m for mat in matches if mat.sids for m in mat.sids]
            synonyms=[m for mat in matches if mat.synonyms for m in mat.synonyms]
            data_frame.loc[index,'PubChem_CID']=';'.join(map(str,cids))
            data_frame.loc[index,'PubChem_SID']=';'.join(map(str,sids[:3]))
            data_frame.loc[index,'PubChem_Synonyms']=';'.join(map(str,synonyms))

        except Exception:
            pass

    return data_frame


class ClassifyMolecules:
    
    def __init__(self,data_frame:pd.DataFrame,col_smiles:str):
        self.data_frame=copy.deepcopy(data_frame)
        self.col_smiles=col_smiles
        self.URL = 'http://classyfire.wishartlab.com'
    
    def _submit_query (self,label:str='Metabolites',query_type='STRUCTURE'):
        
        query_pattern='\n'.join([smi for smi in self.data_frame[self.col_smiles]])
        
        try:
            q = requests.post(f'{self.URL}/queries', json={'label':label, 'query_input': query_pattern, 'query_type':query_type},
                              headers={'Accept': 'application/json', 'Content-Type': 'application/json'})
            
            return q.json()['id']
        except requests.exceptions.HTTPError as e:
            return e.response

    def _get_query(self,query_id:str, data_format:str='json'):
        try:
            r = requests.get(f'{self.URL}/queries/{query_id}.json', headers={'Accept': 'application/json'})
            
        except requests.exceptions.HTTPError as e:
            return e.response
        
        return r.json()
    
    def _add_classification_to_df(self,data_frame,json_data):
        for index in data_frame.index:
            for k, val in json_data['entities'][index].items():
                if 'name' in val:
                    data_frame.loc[index,k]=val['name']
                if 'molecular_framework' in k:
                    data_frame.loc[index,k]=val
                if 'substituents' in k:
                    data_frame.loc[index,k]=';'.join(val)
    
    
    def classify_molecules(self):
        self.query_id=self._submit_query()
        
        self.classification_data=self._get_query(self.query_id)
        
        self._add_classification_to_df(self.data_frame,self.classification_data)
    
    
    
    
    
    
    
    
    '''
    def _get_chemont_node(self,chemontid:str):
        chemont_id = chemontid.replace('CHEMONTID:', 'C')
        try:
            r = requests.get(f'{self.URL}/tax_nodes/{chemont_id}.json', headers={'Accept': 'application/json'})
            return r.json()
        except requests.exceptions.HTTPError as e:
            return e.response
    
    def _get_sequence_classification(self,fingerprint:str, format='json'):
        try:
            if format == 'json':
                r = requests.get(f'{self.URL}/entities/{fingerprint}.{format}', headers={'Accept': 'application/json'})
            else:
                raise ValueError('Invalid format')
        except requests.exceptions.HTTPError as e:
            return e.response
        return r.json()
    '''