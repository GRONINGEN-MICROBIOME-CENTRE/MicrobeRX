"""
This is a module that provides functions to load and process data for the MicrobeRX tool.

The module requires the following packages: pandas, importlib_resources

The module contains the following functions:

- load_reaction_rules: Load the reaction rules from a compressed tab-separated file.
- load_human_evidences: Load the human evidences from a compressed tab-separated file.
- load_microbes_evidences: Load the microbes evidences from a compressed tab-separated file.
- load_microbes_reactions: Load the microbes reactions from a compressed tab-separated file.
- load_microbes_data: Load the microbes data from a compressed tab-separated file.
"""

__all__ = [
    "load_reaction_rules",
    "load_human_evidences",
    "load_microbes_evidences",
    "load_microbes_reactions",
    "load_microbes_data",
]


from importlib_resources import files
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


__REACTION_RULES = files("microberx.DataBase").joinpath("ReactionRules.tsv.gz")
__HUMAN_EVIDENCES = files("microberx.DataBase").joinpath("HumanEvidences.tsv.gz")
__MICROBES_EVIDENCES = files("microberx.DataBase").joinpath("MicrobesEvidences.tsv.gz")
__MICROBES_DATA = files("microberx.DataBase").joinpath("MicrobesData.tsv.gz")
__MICROBES_REACTIONS = files("microberx.DataBase").joinpath("MicrobesReactions.tsv.gz")


def load_reaction_rules():
    """
    Load the reaction rules from a compressed tab-separated file.

    Returns:
        pandas.DataFrame: A dataframe containing the reaction rules, with columns:
            - num_atoms : Number of atoms to match in the query to perfom a prediction. 
            - rule : SMARTS string of the single reactant reaction rule (SRRR).
            - reaction_id : Reaction_id in unified MetaNetX v4.0 id or AGORA2.
            - substrate : MetaNetX id of the Real subtrate of the SRRR.
            - substrate_map : Atom mappeed SMARTS of the of the Real subtrate of the SRRR.  
            - product : MetaNetX id of the Main real subtrate of the SRRR.
            - product_map : Atom mappeed SMARTS of Main real product of the SRRR.
    """
    logging.info("Loading reaction rules...")
    return pd.read_csv(__REACTION_RULES, sep="\t", compression="gzip")


def load_human_evidences():
    """
    Load the human evidences from a compressed tab-separated file.

    Returns:
        pandas.DataFrame: A dataframe containing the human evidences, with columns:
            - source : The unique identifier of the source coming from the metabolic reconstruction.
            - name : Name of the biotransformations, can match with enzyme name.
            - ec : Enzyme Commission number for the biotransformation.
            - mnx_id : Unified id from MetaNetX v4.0. 
            - organisms_count : Number of organims where this souce id has been found.
            - xrefs : coss-references to other reaction databases.
            - origin : Tells if the reaction is coming from human or gut microbes.
            - complexes_count : Numer of genes or complexes found in the metabolic network for this biotransformation. 
    """
    logging.info("Loading human evidences...")
    return pd.read_csv(__HUMAN_EVIDENCES, sep="\t", compression="gzip")


def load_microbes_evidences():
    """
    Load the microbes evidences from a compressed tab-separated file.

    Returns:
        pandas.DataFrame: A dataframe containing the microbes evidences, with columns:
            - source : The unique identifier of the source coming from the metabolic reconstruction.
            - name : Name of the biotransformations, can match with enzyme name.
            - ec : Enzyme Commission number for the biotransformation.
            - mnx_id : Unified id from MetaNetX v4.0. 
            - organisms_count : Number of organims where this souce id has been found.
            - xrefs : coss-references to other reaction databases.
            - origin : Tells if the reaction is coming from human or gut microbes.
            - complexes_count : Numer of genes or complexes found in the metabolic network for this biotransformation. 
    """
    logging.info("Loading microbes evidences...")
    return pd.read_csv(__MICROBES_EVIDENCES, sep="\t", compression="gzip")


def load_microbes_reactions():
    """
    Load the microbes reactions from a compressed tab-separated file.

    Returns:
        pandas.DataFrame: A dataframe containing the microbes reactions.
         - index: strain name of all gut microbes included in microbeRX (source: AGORA2).
         - columns : source name of biotransformation from the metabolic reconstructions.
         - data : any cell contains information about the genes or complexes that have been annotated for each organims and biotransformation.
            
    """
    logging.info("Loading microbes reactions...")
    return pd.read_csv(
        __MICROBES_REACTIONS, sep="\t", index_col=[0], compression="gzip", dtype=str
    )


def load_microbes_data():
    """
    Load the microbes data from a compressed tab-separated file.

    Returns:
        pandas.DataFrame: A dataframe containing the microbes data, with columns:
            - microbe_name
            - Strain
            - Species
            - Genus
            - Family 
            - Order
            - Class
            - Phylum
            - Kingdom
            - Host
            - NCBI Taxonomy ID
            - Cultured
            - Ecosystem
            - Ecosystem Category
            - Ecosystem Subtype
            - Ecosystem Type
            - Gram Staining
            - Oxygen Requirement
            - Motility
    """
    logging.info("Loading microbes data...")
    return pd.read_csv(__MICROBES_DATA, sep="\t", compression="gzip")
