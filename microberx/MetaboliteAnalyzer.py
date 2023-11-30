"""
This module contains functions to compute and analyze molecular properties and descriptors of metabolites using various libraries and web services.

The module has the following functions:

- compute_molecular_descriptors: Computes some molecular descriptors for a given data frame of SMILES strings.
- compute_isotopic_mass: Computes the isotopic mass distribution of a given data frame using the pyOpenMS library.
- search_pubchem: Searches the PubChem database for compounds that match a given data frame of identifiers.
- classify_molecules: Classify molecules based on their SMILES strings using the ClassyFire web service.
"""

__all__ = [
    "compute_molecular_descriptors",
    "compute_isotopic_mass",
    "search_pubchem",
    "classify_molecules"
]


import requests, copy, random, time
import pandas as pd
import numpy as np
import pubchempy as pcp
from pyopenms import *

import plotly.express as px
import plotly.graph_objects as go

from distinctipy import distinctipy

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.FilterCatalog import *

from tqdm.notebook import tqdm


def compute_molecular_descriptors(data_frame: pd.DataFrame, smiles_col: str):
    """ 
    Computes some molecular descriptors and filters for a given data frame of SMILES strings.

    Parameters
    ----------
    data_frame : pd.DataFrame
        A pandas data frame that contains SMILES strings of molecules.
    smiles_col : str
        The name of the column that contains the SMILES strings.

    Returns
    -------
    data_frame : pd.DataFrame
        The same data frame as input, but with additional columns for each molecular descriptor and filter computed. The descriptors and filters are:
            - MolWt: the molecular weight of the molecule
            - LogP: the octanol-water partition coefficient of the molecule
            - NumHAcceptors: the number of hydrogen bond acceptors in the molecule
            - NumHDonors: the number of hydrogen bond donors in the molecule
            - NumRotatableBonds: the number of rotatable bonds in the molecule
            - TPSA: the topological polar surface area of the molecule
            - MolFormula: the molecular formula of the molecule
            - Lipinski: a boolean value that indicates whether the molecule satisfies the Lipinski's rule of five or not. The rule of five states that most drug-like molecules have molecular weight less than 500, LogP less than 5, number of hydrogen bond acceptors less than 10, and number of hydrogen bond donors less than 5.
            - Veber: a boolean value that indicates whether the molecule satisfies the Veber's rule or not. The rule states that most orally active drugs have 10 or fewer rotatable bonds and a polar surface area equal to or less than 140 Å2.
            - Brenk: a string that contains the names of the Brenk filters that the molecule matches, separated by semicolons. The Brenk filters are a set of 68 unwanted substructures that are associated with reactive fucntional groups.
            - PAINS: a string that contains the names of the PAINS filters that the molecule matches, separated by semicolons. The PAINS filters are a set of 480 substructures that are associated with pan-assay interference compounds.
    """
    brenk_params = FilterCatalogParams()
    brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    brenk = FilterCatalog(brenk_params)

    pains_params = FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains = FilterCatalog(pains_params)
    
    data_frame = copy.deepcopy(data_frame)
    
    for index in data_frame.index:
        mol = Chem.MolFromSmiles(data_frame[smiles_col][index])
        if mol:
            data_frame.loc[index, "MolWt"] = round(Descriptors.MolWt(mol), 3)
            data_frame.loc[index, "LogP"] = round(Descriptors.MolLogP(mol), 3)
            data_frame.loc[index, "NumHAcceptors"] = Descriptors.NumHAcceptors(mol)
            data_frame.loc[index, "NumHDonors"] = Descriptors.NumHDonors(mol)
            data_frame.loc[index, "NumRotatableBonds"] = Descriptors.NumRotatableBonds(
                mol
            )
            data_frame.loc[index, "TPSA"] = Descriptors.TPSA(mol)
            data_frame.loc[index, "MolFormula"] = AllChem.CalcMolFormula(mol)

            if data_frame["NumHDonors"][index] > 5 or data_frame["NumHAcceptors"][index] > 10 or data_frame["MolWt"][index] > 500 or data_frame["LogP"][index] > 5 :
                data_frame.loc[index,"Lipinski"]=False
            else: 
                data_frame.loc[index,"Lipinski"]=True
        
            if data_frame["NumHDonors"][index] > 5 or data_frame["NumHAcceptors"][index] > 10 or data_frame["MolWt"][index] > 500 or data_frame["LogP"][index] > 5 or data_frame["NumRotatableBonds"][index] > 10 or data_frame["TPSA"][index] > 140:
                data_frame.loc[index,"Veber"]=False
            else:
                data_frame.loc[index,"Veber"]=True
        
            brenk_results=brenk.GetMatches(mol)
            pains_results=pains.GetMatches(mol)

            if brenk_results:
                data_frame.loc[index,"Brenk"]=";".join([r.GetDescription() for r in brenk_results])
            if pains_results:
                data_frame.loc[index,"PAINS"]=";".join([r.GetDescription() for r in pains_results])
                
    return data_frame


def compute_isotopic_mass(data_frame: pd.DataFrame, molformula_col: str):
    """
    Computes the isotopic mass distribution of a given data frame using the pyOpenMS library.

    Parameters
    ----------
    data_frame : pd.DataFrame
        A pandas data frame that contains the molecular formulas as a column.
    molformula_col : str
        A string that specifies the name of the column that contains the molecular formulas.

    Returns
    -------
    data_frame : pd.DataFrame
        A pandas data frame that has two additional columns: 'probability_sum' and 'mass_distribution'. The 'probability_sum' column contains the sum of the probabilities of all isotopes for each molecular formula. The 'mass_distribution' column contains the mass and probability of each isotope as a string, separated by semicolons.

    Example
    -------
    >>> import pandas as pd
    >>> from pyopenms import EmpiricalFormula, CoarseIsotopePatternGenerator
    >>> df = pd.DataFrame({'formula': ['C6H12O6', 'C2H4O2', 'C3H8O3']})
    >>> df = compute_isotopic_mass(df, 'formula')
    >>> print(df)
        formula  probability_sum                                   mass_distribution
    0  C6H12O6            1.0000  180.0634:100.0;181.0668:10.72;182.0701:1.176;183...
    1   C2H4O2            1.0000  60.0211:100.0;61.0245:11.08;62.0279:1.216;63.031...
    2   C3H8O3            0.9999  92.0473:100.0;93.0507:10.55;94.0541:1.159;95.057...
    """
    data_frame = copy.deepcopy(data_frame)
    for index in data_frame.index:
        try:
            molF = EmpiricalFormula(data_frame[molformula_col][index])
            isotopes = molF.getIsotopeDistribution(CoarseIsotopePatternGenerator(4))
            prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
            distribution = ";".join(
                [
                    f"{round(iso.getMZ(),4)}:{round(iso.getIntensity()*100,4)}"
                    for iso in isotopes.getContainer()
                ]
            )
            sumR = round(prob_sum, 4)
            data_frame.loc[index, "probability_sum"] = sumR
            data_frame.loc[index, "mass_distribution"] = distribution
        except Exception:
            pass

    return data_frame


def search_pubchem(data_frame: pd.DataFrame, entry_col: str, entry_type: str = "smiles"):
    """
    Searches the PubChem database for compounds that match a given data frame of identifiers.

    Parameters
    ----------
    data_frame : pd.DataFrame
        A pandas data frame that contains the identifiers of the compounds to search for.
    entry_col : str
        A string that specifies the name of the column that contains the identifiers.
    entry_type : str, optional
        A string that specifies the type of the identifiers, such as 'smiles', 'inchi', 'cid', etc. The default is 'smiles'.

    Returns
    -------
    data_frame : pd.DataFrame
        The same data frame as input, but with additional columns for each PubChem property retrieved. The properties are:
            - PubChem_CID: the PubChem compound identifier, separated by semicolons if there are multiple matches.
            - PubChem_SID: the PubChem substance identifier, separated by semicolons if there are multiple matches. Only the first three SIDs are shown.
            - PubChem_Synonyms: the synonyms of the compound, separated by semicolons if there are multiple matches.
    """
    for index in tqdm(data_frame.index):
        try:
            matches = pcp.get_compounds(
                data_frame[entry_col][index], namespace=entry_type
            )
            cids = [mat.cid for mat in matches]
            sids = [m for mat in matches if mat.sids for m in mat.sids]
            synonyms = [m for mat in matches if mat.synonyms for m in mat.synonyms]
            data_frame.loc[index, "PubChem_CID"] = ";".join(map(str, cids))
            data_frame.loc[index, "PubChem_SID"] = ";".join(map(str, sids[:3]))
            data_frame.loc[index, "PubChem_Synonyms"] = ";".join(map(str, synonyms))

        except Exception:
            pass

    return data_frame


def classify_molecules(data_frame: pd.DataFrame, smiles_col: str, names_col: str):
    """
    Classify molecules based on their SMILES strings.

    This function submits a query to the ClassyFire web service and returns a data frame with the classification results.

    Parameters
    ----------
    data_frame : pd.DataFrame
        The input data frame with the molecules information.
    smiles_col : str
        The name of the column that contains the SMILES strings.
    names_col : str
        The name of the column that contains the molecule names.

    Returns
    -------
    pd.DataFrame
        The output data frame with the classification results added as new columns. The columns are:
            - kingdom: the name of the chemical kingdom of the molecule, such as 'Organic compounds', 'Inorganic compounds', etc.
            - superclass: the name of the chemical superclass of the molecule, such as 'Lipids and lipid-like molecules', 'Organoheterocyclic compounds', etc.
            - class: the name of the chemical class of the molecule, such as 'Steroids and steroid derivatives', 'Benzodiazepines', etc.
            - subclass: the name of the chemical subclass of the molecule, such as 'Cholestane steroids', '1,4-benzodiazepines', etc.

    Raises
    ------
    requests.exceptions.HTTPError
        If the query to the ClassyFire web service fails.
    """
    
    URL = "http://classyfire.wishartlab.com"

    def _submit_query(data_frame: pd.DataFrame,smiles_col: str,names_col: str,label: str = "MicrobeRX",query_type:str="STRUCTURE"):
        unique_mols = data_frame.drop_duplicates(subset=[names_col, smiles_col])

        entries = [
            f"{unique_mols[names_col][index]}\t{unique_mols[smiles_col][index]}"
            for index in unique_mols.index
        ]

        query_pattern = "\n".join(entries)
        
        try:
            q = requests.post(
                f"{URL}/queries",
                json={
                    "label": label,
                    "query_input": query_pattern,
                    "query_type": query_type,
                    "fstruc_content_type": "tsv",
                },
                headers={
                    "Accept": "application/json",
                    "Content-Type": "application/json",
                },
            )
            
            return q.json()["id"]
        except requests.exceptions.HTTPError as e:
            return e.response


    def _get_query(query_id: str, data_format: str = "json"):
        try:
            r = requests.get(
                f"{URL}/queries/{query_id}.json", headers={"Accept": "application/json"}
            )

        except requests.exceptions.HTTPError as e:
            return e.response
    
        return r.json()

    
    data_frame = copy.deepcopy(data_frame)
    
    if len(data_frame.index)>100:
        chunks=np.array_split(data_frame, round(len(data_frame.index)/100))
    else:
        chunks=[data_frame]
    
    job_ids=[]
    for c in chunks:
        job=_submit_query(c,smiles_col=smiles_col,names_col=names_col)
        job_ids.append(job)
    
    time.sleep(10)
    
    results=[]
    for result in job_ids:
        data=_get_query(str(result))
        results.append(pd.DataFrame(data['entities']))
    
    class_data=pd.concat(results,ignore_index=True)
    
    class_data=class_data.transform(lambda x: x.apply(lambda y: y['name'] if isinstance(y,dict) else y))
    class_data['alternative_parents']=class_data['alternative_parents'].apply(lambda x: ';'.join([y['name'] for y in x] if isinstance(x,list) else np.nan))
    class_data['intermediate_nodes']=class_data['intermediate_nodes'].apply(lambda x: ';'.join([y['name'] for y in x] if isinstance(x,list) else np.nan))
    class_data['external_descriptors']=class_data['external_descriptors'].apply(lambda x: ';'.join([f"{y['source']}:{y['source_id']}" for y in x] if isinstance(x,list) else np.nan))
    class_data=class_data.transform(lambda x: x.apply(lambda y: ';'.join(y) if isinstance(y,list) else y))

    data_frame_classification=data_frame.merge(class_data,left_on=names_col,right_on='identifier', how='outer')
    
    return data_frame_classification