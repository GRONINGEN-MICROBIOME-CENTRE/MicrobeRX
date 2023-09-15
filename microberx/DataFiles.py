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

__version__ = "0.1.4"

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
            - rule_id: The unique identifier of the rule
            - reactant: The chemical formula of the reactant
            - product: The chemical formula of the product
            - rate_constant: The rate constant of the reaction
    """
    logging.info("Loading reaction rules...")
    return pd.read_csv(__REACTION_RULES, sep="\t", compression="gzip")


def load_human_evidences():
    """
    Load the human evidences from a compressed tab-separated file.

    Returns:
        pandas.DataFrame: A dataframe containing the human evidences, with columns:
            - evidence_id: The unique identifier of the evidence
            - gene: The gene symbol of the gene involved in the evidence
            - disease: The disease name of the disease involved in the evidence
            - score: The score of the evidence, ranging from 0 to 1
    """
    return pd.read_csv(__HUMAN_EVIDENCES, sep="\t", compression="gzip")


def load_microbes_evidences():
    """
    Load the microbes evidences from a compressed tab-separated file.

    Returns:
        pandas.DataFrame: A dataframe containing the microbes evidences, with columns:
            - evidence_id: The unique identifier of the evidence
            - microbe: The scientific name of the microbe involved in the evidence
            - disease: The disease name of the disease involved in the evidence
            - score: The score of the evidence, ranging from 0 to 1
    """
    logging.info("Loading human evidences...")
    return pd.read_csv(__MICROBES_EVIDENCES, sep="\t", compression="gzip")


def load_microbes_reactions():
    """
    Load the microbes reactions from a compressed tab-separated file.

    Returns:
        pandas.DataFrame: A dataframe containing the microbes reactions, with columns:
            - reaction_id: The unique identifier of the reaction
            - microbe: The scientific name of the microbe involved in the reaction
            - reactant: The chemical formula of the reactant
            - product: The chemical formula of the product
            - rate_constant: The rate constant of the reaction
    """
    logging.info("Loading microbes evidences...")
    return pd.read_csv(
        __MICROBES_REACTIONS, sep="\t", index_col=[0], compression="gzip", dtype=str
    )


def load_microbes_data():
    """
    Load the microbes data from a compressed tab-separated file.

    Returns:
        pandas.DataFrame: A dataframe containing the microbes data, with columns:
            - microbe_id: The unique identifier of the microbe
            - name: The scientific name of the microbe
            - type: The type of the microbe (e.g. bacteria, virus, fungus, etc.)
            - genome: The genome sequence of the microbe
    """
    logging.info("Loading microbes reactions...")
    return pd.read_csv(__MICROBES_DATA, sep="\t", compression="gzip")
