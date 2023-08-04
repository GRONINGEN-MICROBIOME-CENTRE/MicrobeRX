"""
Location of data files for the MicrobeRX databases
====================================================

"""

from importlib_resources import files

__all__ = ["load_reaction_rules", "load_human_evidences", 
           "load_microbes_evidences", "load_microbes_reactions", "load_microbes_data"]

__REACTION_RULES      = files("MicrobeRX.DataBase").joinpath("ReactionRules.tsv.gz")
__HUMAN_EVIDENCES     = files("MicrobeRX.DataBase").joinpath("HumanEvidences.tsv.gz")
__MICROBES_EVIDENCES  = files("MicrobeRX.DataBase").joinpath("MicrobesEvidences.tsv.gz")
__MICROBES_DATA       = files("MicrobeRX.DataBase").joinpath("MicrobesData.tsv.gz")
__MICROBES_REACTIONS  = files("MicrobeRX.DataBase").joinpath("MicrobesReactions.tsv.gz")

''' 
################### DATABASES ###################
'''
import pandas as pd

def load_reaction_rules():    
    return pd.read_csv(__REACTION_RULES,sep='\t',compression='gzip')

def load_human_evidences():
    return pd.read_csv(__HUMAN_EVIDENCES,sep='\t',compression='gzip')

def load_microbes_evidences():
    return pd.read_csv(__MICROBES_EVIDENCES,sep='\t',compression='gzip')

def load_microbes_reactions():    
    return pd.read_csv(__MICROBES_REACTIONS,sep='\t',index_col=[0],compression='gzip')

def load_microbes_data():    
    return pd.read_csv(__MICROBES_DATA,sep='\t',compression='gzip')

