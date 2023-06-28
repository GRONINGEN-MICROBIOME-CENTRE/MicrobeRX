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


class Analyzer:
    
    def __init__(self,data_frame:pd.DataFrame,col_smiles:str):
        self.data_frame=copy.deepcopy(data_frame)
        self.col_smiles=col_smiles
    
    def compute_molecular_descriptors(self):
        
        self.molecules=[]
        for index in tqdm(self.data_frame.index):
            mol=Chem.MolFromSmiles(self.data_frame[self.col_smiles][index])
            if mol:
                self.molecules.append(mol)

                self.data_frame.loc[index,'MolWt']= round(Descriptors.MolWt(mol),3) 
                self.data_frame.loc[index,'LogP']= round(Descriptors.MolLogP(mol),3) 
                self.data_frame.loc[index,'NumHAcceptors']= Descriptors.NumHAcceptors(mol) 
                self.data_frame.loc[index,'NumHDonors']= Descriptors.NumHDonors(mol) 
                self.data_frame.loc[index,'NumRotatableBonds']= Descriptors.NumRotatableBonds(mol)
                self.data_frame.loc[index,'TPSA']= Descriptors.TPSA(mol) 
                self.data_frame.loc[index,'MolFormula']= AllChem.CalcMolFormula(mol)
                
    def plot_molecular_descriptors(self) -> go.Figure:
        """
        Plot molecular descriptors values in a polar plot.

        Parameters:
        - table: pd.DataFrame: Dataframe containing the molecular descriptors.
        - outfile: str: File name to save the figure to, as an image. If None, the function will return the figure object.

        Returns:
        - F: go.Figure: Figure object of the plotly polar plot.

        """
        data=pd.DataFrame(index=list(self.data_frame.metabolite_id)) 
        data['MolWt']=[i/500 for i in self.data_frame['MolWt']]
        data['LogP']=[i/5 for i in self.data_frame['LogP']]
        data['HBA']=[i/10 for i in self.data_frame['NumHAcceptors']]
        data['HBD']=[i/5 for i in self.data_frame['NumHDonors']]
        data['RotB']=[i/10 for i in self.data_frame['NumRotatableBonds']]
        data['TPSA']=[i/140 for i in self.data_frame['TPSA']]

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
        
    def compute_isotopic_mass(self):
          
        for index in tqdm(self.data_frame.index):
            try:
                molF = EmpiricalFormula(self.data_frame['MolFormula'][index])
                isotopes=molF.getIsotopeDistribution(CoarseIsotopePatternGenerator(4))
                prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
                distribution=';'.join([f'{round(iso.getMZ(),4)}:{round(iso.getIntensity()*100,4)}' for iso in isotopes.getContainer()])
                sumR=round(prob_sum,4)  
            except Exception:
                sumR=None
                distribution=None
                
            self.data_frame.loc[index,'probability_sum']=sumR
            self.data_frame.loc[index,'mass_distribution']=distribution
            index=index+1
            
    def plot_isotopic_masses(self):
        
        masses=copy.deepcopy(self.data_frame[['main_product_smiles','metabolite_id','mass_distribution']])
        masses.dropna(subset=['mass_distribution'],inplace=True)
        masses.mass_distribution=masses.mass_distribution.apply(lambda x: {float(value.split(':')[0]):float(value.split(':')[1]) for value in x.split(';')})
        masses.reset_index(drop=True,inplace=True)
        
        Figure = go.Figure()

        color = [distinctipy.get_hex(c) for c in distinctipy.get_colors(n_colors=len(masses.index),pastel_factor=0.7)]
        
        for index in masses.index:
            f=px.bar(x=tuple(masses.mass_distribution[index].keys()), y=tuple(masses.mass_distribution[index].values()),color_discrete_sequence=[color[index]])
            f.update_traces(name=masses.metabolite_id[index])
            f.data[0].showlegend=True
            Figure.add_trace(f.data[0])

        Figure.update_layout(width=1200,height=700,showlegend=True,legend=dict(x=1.2, y=0.95),legend_orientation="h",plot_bgcolor='rgba(0, 0, 0, 0)',font=dict(size=12),xaxis_title="Atomic Mass (u)",yaxis_title="Relative abundance (%)",title="Isotopic distribution")

        return Figure
    
    def search_pubchem(self,entry_type:str='smiles'):
        
        for index in tqdm(self.data_frame.index):
            try:
                matches=pcp.get_compounds(self.data_frame[self.col_smiles][index],namespace=entry_type)
                cids=[mat.cid for mat in matches]
                sids=[m for mat in matches if mat.sids for m in mat.sids]
                synonyms=[m for mat in matches if mat.synonyms for m in mat.synonyms]
                self.data_frame.loc[index,'PubChem_CID']=';'.join(map(str,cids))
                self.data_frame.loc[index,'PubChem_SID']=';'.join(map(str,sids[:3]))
                self.data_frame.loc[index,'PubChem_Synonyms']=';'.join(map(str,synonyms))

            except Exception:
                pass
            
            index=index+1
        
    

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