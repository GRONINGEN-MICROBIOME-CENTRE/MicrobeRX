"""
Location of data files for the MicrobeRX unit tests
====================================================

"""

from importlib_resources import files

__all__ = ["REACTION_RULES", "HUMAN_EVIDENCES", "MICROBES_EVIDENCES", "MICROBES_DATA"]

REACTION_RULES      = files("package.MicrobeRX.DataBase").joinpath("ReactionRules.tsv.gz")
HUMAN_EVIDENCES     = files("package.MicrobeRX.DataBase").joinpath("HumanEvidences.tsv.gz")
MICROBES_EVIDENCES  = files("package.MicrobeRX.DataBase").joinpath("MicrobesEvidences.tsv.gz")
MICROBES_DATA       = files("package.MicrobeRX.DataBase").joinpath("MicrobesData.tsv.gz")