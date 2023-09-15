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


def compute_molecular_descriptors(data_frame: pd.DataFrame, smiles_col: str):
    """
    Computes some molecular descriptors for a given data frame of SMILES strings.

    Parameters
    ----------
    data_frame : pd.DataFrame
        A pandas data frame that contains SMILES strings of molecules.
    smiles_col : str
        The name of the column that contains the SMILES strings.

    Returns
    -------
    data_frame : pd.DataFrame
        The same data frame as input, but with additional columns for each molecular descriptor computed. The descriptors are:
            - MolWt: the molecular weight of the molecule
            - LogP: the octanol-water partition coefficient of the molecule
            - NumHAcceptors: the number of hydrogen bond acceptors in the molecule
            - NumHDonors: the number of hydrogen bond donors in the molecule
            - NumRotatableBonds: the number of rotatable bonds in the molecule
            - TPSA: the topological polar surface area of the molecule
            - MolFormula: the molecular formula of the molecule
    """
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

    def _submit_query(data_frame: pd.DataFrame,smiles_col: str,names_col: str,label: str = "Metabolites",query_type="STRUCTURE"):
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
                    "label": "caca",
                    "query_input": query_pattern,
                    "query_type": "STRUCTURE",
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

    def _add_classification_to_df(data_frame, json_data, names_col):
        data = {}
        for e in json_data["entities"]:
            cl = {}
            for k, v in e.items():
                if "name" in v:
                    cl[k] = v["name"]
                if "molecular_framework" in k:
                    cl[k] = v
                if "substituents" in k:
                    cl[k] = ";".join(v)
            data[e["identifier"]] = cl

        for index in data_frame.index:
            try:
                indx_class = data[data_frame[names_col][index]]
                data_frame.loc[index, indx_class.keys()] = indx_class.values()
            except Exception:
                pass

        return data_frame

    data_frame = copy.deepcopy(data_frame)

    job_id = _submit_query(data_frame, smiles_col, names_col)
    results = _get_query(job_id)
    data_frame_classification = _add_classification_to_df(
        data_frame, results, names_col
    )

    return data_frame_classification
