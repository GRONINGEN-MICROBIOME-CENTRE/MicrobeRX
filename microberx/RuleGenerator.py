"""
This is a module that provides functions to manipulate chemical reactions and and generate reaction rules.

The module contains the following functions:

- decompose_reaction: A fucntion to parse chemical reaction from a string format into a dictionary format.

- sanitize_reaction: Sanitizes a chemical reaction by removing stereochemistry, replacing dummy atoms with carbon, and standardizing the molecules.

- map_reaction: Maps the atoms of a chemical reaction using either RXNMapper or ReactionDecoder.

- set_reaction_ids: Sets the IDs of the reactants and products of a target reaction based on their molecular formulas and a reference reaction.

- reverse_reaction: Reverses a chemical reaction by swapping the reactants and products.

- generate_single_reactant_reactions: Generates a dictionary of single reactant reactions from a mapped reaction.

- generate_rules: Generates a dictionary of rules for a single reactant reaction based on the reacting atoms and the rings.

- Reaction: A class to represent a chemical reaction with various attributes and methods.

"""


__all__ = [
    "decompose_reaction",
    "sanitize_reaction",
    "sanitize_reaction",
    "set_reaction_ids",
    "reverse_reaction",
    "generate_single_reactant_reactions",
    "generate_rules",
    "Reaction",
]


import copy, subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, rdqueries, DataStructs, MACCSkeys

import datamol as dm

from rxnmapper import RXNMapper

rxn_mapper = RXNMapper()

import pandas as pd
import re

from importlib_resources import files

REACTION_DECODER = files("microberx.bin").joinpath("RTD.jar")

def decompose_reaction(reaction: str = None, compounds_map: dict = None):
    """
    Decompose a chemical reaction into its individual compounds and stoichiometry.

    Parameters
    ----------
    reaction : str, optional
        The chemical reaction in string format. It can be specified using various notations:
        - "<=>" for reversible reactions.
        - "-->" for irreversible reactions.
        - "=" for generic reactions.

    compounds_map : dict, optional
        A dictionary mapping compound names to their corresponding SMILES representations.

    Returns
    -------
    reaction_dict : dict
        A dictionary containing the decomposed information of the chemical reaction, including:
        - "Reaction": The original input chemical reaction.
        - "Reversible": A boolean indicating whether the reaction is reversible.
        - "LEFT": A dictionary of reactants with stoichiometry and SMILES.
        - "RIGHT": A dictionary of products with stoichiometry and SMILES.
        - "ReactionSmiles": The SMILES representation of the entire reaction.
        - "ReactionNames": The names of compounds in the reaction.

    Examples
    --------
    >>> compounds_map = {"H2O": "O", "CO2": "O=C=O", "CH4": "C"}
    >>> reaction_str = "H2O + CO2 <=> CH4 + H2O"
    >>> reaction_info = decompose_reaction(reaction_str, compounds_map)
    >>> print(reaction_info)
    >>> reaction_str = "H2O + CO2 --> CH4 + H2O"
    >>> reaction_info = decompose_reaction(reaction_str, compounds_map)
    >>> print(reaction_info)
    """
    
    if " = " in reaction:
        reaction = re.sub("@\w+", "", reaction)
        reactants, products = reaction.split(" = ")
        Reversible = True

    if "<=>" in reaction:
        reaction = re.sub(r"\[[a-z]\]", "", reaction)
        reactants, products = reaction.split("<=>")
        Reversible = True

    if "-->" in reaction:
        reaction = re.sub(r"\[[a-z]\]", "", reaction)
        reactants, products = reaction.split("-->")
        Reversible = False

    # Split the reactants and products into individual compounds
    reactants = reactants.split("+")
    products = products.split("+")

    # Create empty dictionaries to store the reactants and products with their stoichiometry
    Left_side = {}
    Right_side = {}
    smiles_scheme = ""
    names_scheme = ""

    for r in reactants:
        r = r.strip()
        if len(r.split()) > 1:
            stoich, name = r.split()
        else:
            name = r
            stoich = 1

        smiles = compounds_map[name]
        Left_side[name] = {"stoichiometry": float(stoich), "smiles": smiles}
        smiles_scheme += (smiles + ".") * int(float(stoich))
        names_scheme += (name + ".") * int(float(stoich))

    smiles_scheme = smiles_scheme[:-1]
    names_scheme = names_scheme[:-1]

    smiles_scheme += ">>"
    names_scheme += ">>"
    for p in products:
        p = p.strip()
        if len(p.split()) > 1:
            stoich, name = p.split()
        else:
            name = p
            stoich = 1

        smiles = compounds_map[name]
        Right_side[name] = {"stoichiometry": float(stoich), "smiles": smiles}
        smiles_scheme += (smiles + ".") * int(float(stoich))
        names_scheme += (name + ".") * int(float(stoich))

    reaction_dict = {
        "Reaction": reaction,
        "Reversible": Reversible,
        "LEFT": Left_side,
        "RIGHT": Right_side,
        "ReactionSmiles": smiles_scheme[:-1],
        "ReactionNames": names_scheme[:-1],
    }

    return reaction_dict


def sanitize_reaction(target_reaction: AllChem.ChemicalReaction):
    """
    Sanitizes a chemical reaction by removing stereochemistry, replacing dummy atoms with carbon, and standardizing the molecules.

    Parameters
    ----------
    target_reaction : AllChem.ChemicalReaction
        The input chemical reaction to be sanitized.

    Returns
    -------
    fixed_reaction : AllChem.ChemicalReaction
        The output chemical reaction after sanitization. The fixed reaction has the following features:
            - The stereochemistry of the reactants and products is removed, as it may not be relevant or accurate for the reaction.
            - The dummy atoms (*) in the reactants and products are replaced with carbon atoms (#6), as they may represent unspecified groups or atoms.
            - The reactants and products are sanitized and standardized using the dm module, which performs operations such as kekulization, neutralization, tautomerization, etc.
            - The 2D coordinates of the reactants and products are computed using the AllChem.Compute2DCoords function, which may improve the visualization of the reaction.
    Examples
    --------
    >>>reaction = AllChem.ReactionFromSmarts("[OH:1].[C:2]=O>>[C:2][OH:1]")
    >>>fixed_reaction = sanitize_reaction(reaction)
    >>>img = Draw.ReactionToImage(fixed_reaction)
    >>>img.show()
    """

    def _replace_dummy_atoms(mol: Chem.Mol):
        """
        Performs a molecular sanitization, removes stereochemistry, and replaces dummy atoms.

        Parameters
        ----------
        mol : Chem.Mol
            The input molecule to be sanitized.

        Returns
        -------
        fixed_mol : Chem.Mol
            The output molecule after sanitization.
        """
        dummy_query = Chem.MolFromSmiles("*")
        fixed_mol = AllChem.ReplaceSubstructs(
            mol,
            dummy_query,
            Chem.MolFromSmiles("[#6]"),
            replaceAll=True,
            useChirality=False,
        )[0]
        fixed_mol = dm.fix_mol(fixed_mol)
        fixed_mol = dm.sanitize_mol(fixed_mol)
        fixed_mol = dm.standardize_mol(fixed_mol)
        Chem.RemoveStereochemistry(fixed_mol)
        AllChem.Compute2DCoords(fixed_mol)
        return fixed_mol

    fixed_reaction = AllChem.ChemicalReaction()

    for index, reactant in enumerate(target_reaction.GetReactants()):
        fixed_reactant = _replace_dummy_atoms(reactant)
        if fixed_reactant != None:
            fixed_reaction.AddReactantTemplate(fixed_reactant)
    for index, product in enumerate(target_reaction.GetProducts()):
        fixed_product = _replace_dummy_atoms(product)
        if fixed_product != None:
            fixed_reaction.AddProductTemplate(fixed_product)

    return fixed_reaction


def map_reaction(reaction_smiles: str = None, mapper: str = "RXNMapper"):
    """
    Maps the atoms of a chemical reaction using either RXNMapper or ReactionDecoder.

    Parameters
    ----------
    reaction_smiles : str, optional
        The input chemical reaction in SMILES format. Default is None.
    mapper : str, optional
        The name of the mapper to use. Either 'RXNMapper' or 'ReactionDecoder'. Default is 'RXNMapper'.

    Returns
    -------
    mapped_reaction : AllChem.ChemicalReaction
        The output rdkit chemical reaction with atom mapping. The mapped reaction has the following features:
            - The reactants and products are the same as the input reaction, but with atom numbers added to each atom symbol.
            - The atom numbers indicate the correspondence between the atoms in the reactants and products, according to the mapper algorithm.
            - The atom numbers are enclosed in square brackets and separated by colons, such as [1:2], meaning that atom 1 in the reactant becomes atom 2 in the product.
            - The atom numbers are consistent with the SMILES notation, meaning that they follow the order of appearance of the atoms in the SMILES string.
            - The bond types and stereochemistry of the reaction are preserved in the mapped reaction.
    """

    def _reaction_decoder(reaction_smiles: str):
        """
        Internal function to map a reaction using ReactionDecoder.

        Parameters
        ----------
        reaction_smiles : str
            The input chemical reaction in SMILES format.

        Returns
        -------
        mapped_rxn : AllChem.ChemicalReaction
            The mapped reaction using ReactionDecoder.
        """
        out = subprocess.call(
            ["java",
             "-jar",
             REACTION_DECODER,
             "-Q",
             "SMI",
             "-q",
             reaction_smiles,
             "-c",
             "-j",
             "AAM",
             "-f",
             "TEXT"],
            timeout=360)
        
        mapped_rxn = AllChem.ReactionFromRxnFile("ECBLAST_smiles_AAM.rxn")
        return mapped_rxn

    if mapper == "RXNMapper":
        mapper = rxn_mapper.get_attention_guided_atom_maps([reaction_smiles])[0]

        mapped_reaction = AllChem.ReactionFromSmarts(
            mapper["mapped_rxn"], useSmiles=True
        )

    if mapper == "ReactionDecoder":
        mapped_reaction = _reaction_decoder(reaction_smiles)

    return mapped_reaction


def set_reaction_ids(reference_reaction: AllChem.ChemicalReaction,target_reaction: AllChem.ChemicalReaction,reaction_ids: str):
    
    """
    Sets the IDs of the reactants and products of a target reaction based on their molecular formulas and a reference reaction.

    Parameters
    ----------
    reference_reaction : AllChem.ChemicalReaction
        The reference chemical reaction that has the same reactants and products as the target reaction, but in a different order or orientation.
    target_reaction : AllChem.ChemicalReaction
        The target chemical reaction that needs to have its IDs set.
    reaction_ids : str
        The IDs of the reactants and products of the reference reaction, in the format of 'R1.R2.>>P1.P2.', where R1, R2, P1, P2 are the IDs.

    Returns
    -------
    target_reaction : AllChem.ChemicalReaction
        The target chemical reaction with its IDs set according to the reference reaction and the reaction_ids. The target reaction has the following features:
            - The reactants and products are the same as the input target reaction, but with an additional property 'ID' added to each molecule.
            - The value of the 'ID' property is determined by matching the molecular formula of each molecule in the target reaction with the corresponding molecule in the reference reaction, and then using the value from the reaction_ids string.
            - The order and orientation of the reactants and products in the target reaction are preserved in the output target reaction.

    Examples
    --------
    >>> reference_reaction = Chem.ReactionFromSmarts("[H][C:1]([H])=[O:2]>>[O:2]=[C:1]([H])")
    >>> target_reaction = Chem.ReactionFromSmarts("[O:2]=[C:1]([H])>>[H][C:1]([H])=[O:2]")
    >>> reaction_ids = "R1.R2.>>P1.P2."
    >>> target_reaction_with_ids = set_reaction_ids(reference_reaction, target_reaction, reaction_ids)
    >>> print(target_reaction_with_ids.GetReactants()[0].GetProp("ID"))
    >>> print(target_reaction_with_ids.GetProducts()[0].GetProp("ID"))
    """
    
    reactants_ids = reaction_ids.split(">>")[0].split(".")
    products_ids = reaction_ids.split(">>")[1].split(".")

    reactants_formulas = [
        AllChem.CalcMolFormula(r) for r in reference_reaction.GetReactants()
    ]
    products_formulas = [
        AllChem.CalcMolFormula(r) for r in reference_reaction.GetProducts()
    ]

    reactants_map = dict(zip(reactants_formulas, reactants_ids))
    products_map = dict(zip(products_formulas, products_ids))

    for r in target_reaction.GetReactants():
        r.SetProp("ID", reactants_map[AllChem.CalcMolFormula(r)])

    for p in target_reaction.GetProducts():
        p.SetProp("ID", products_map[AllChem.CalcMolFormula(p)])

    return target_reaction


def reverse_reaction(reaction: AllChem.ChemicalReaction):
    """
    Reverses a chemical reaction by swapping the reactants and products.

    Parameters
    ----------
    reaction : AllChem.ChemicalReaction
        The input chemical reaction to be reversed.

    Returns
    -------
    reversed_reaction : AllChem.ChemicalReaction
        The output chemical reaction that is the reverse of the input reaction. The reversed reaction has the following features:
            - The reactants are the same as the products of the input reaction, but in the opposite order.
            - The products are the same as the reactants of the input reaction, but in the opposite order.
            - The atom mapping, bond types, and stereochemistry of the reaction are preserved in the reversed reaction.

    Examples
    --------
    >>> reaction = Chem.ReactionFromSmarts("[H][C:1]([H])=[O:2]>>[O:2]=[C:1]([H])")
    >>> reversed_reaction = reverse_reaction(reaction)
    >>> print(Chem.MolToSmarts(reversed_reaction))
    """

    reversed_reaction = AllChem.ChemicalReaction()

    for i, p in enumerate(reaction.GetProducts()):
        reversed_reaction.AddReactantTemplate(p)

    for i, r in enumerate(reaction.GetReactants()):
        reversed_reaction.AddProductTemplate(r)

    return reversed_reaction


def generate_single_reactant_reactions(mapped_reaction: AllChem.ChemicalReaction):
    """
    Generates a dictionary of single reactant reactions from a mapped reaction.

    Parameters
    ----------
    mapped_reaction:AllChem.ChemicalReaction
        The input chemical reaction with atom mapping.

    Returns
    -------
    all_unique_reactions: dict
        A dictionary of single reactant reactions, keyed by the reactant index and containing the reactant ID and the single reactant reaction object.
    
    Examples
    --------
    >>> mapped_reaction = Chem.ReactionFromSmarts("[H][C:1]([H])=[O:2]>>[O:2]=[C:1]([H])")
    >>> single_reactant_reactions = generate_single_reactant_reactions(mapped_reaction)
    >>> print(single_reactant_reactions["reactantIdx_1"]["ID"])
    >>> print(Chem.MolToSmarts(single_reactant_reactions["reactantIdx_1"]["SingleReactantReaction"]))
    """

    def _sort_reaction(
        single_reactant_reaction: AllChem.ChemicalReaction,
    ) -> AllChem.ChemicalReaction:
        reactant_atom_map = [
            atom.GetAtomMapNum()
            for atom in single_reactant_reaction.GetReactantTemplate(0).GetAtoms()
        ]

        if len(reactant_atom_map) > 0:
            sorted_reaction = AllChem.ChemicalReaction()
            [
                sorted_reaction.AddReactantTemplate(reactant)
                for reactant in single_reactant_reaction.GetReactants()
            ]
            products = {}
            for product in single_reactant_reaction.GetProducts():
                product_atom_map = [
                    atom.GetAtomMapNum()
                    for atom in product.GetAtoms()
                    if atom.GetAtomMapNum() != 0
                ]
                products[product] = len(set(product_atom_map) & set(reactant_atom_map))

            sorted_products = dict(
                sorted(products.items(), key=lambda x: x[1], reverse=True)
            )

            [
                sorted_reaction.AddProductTemplate(product)
                for product in sorted_products.keys()
            ]

        return sorted_reaction

    all_unique_reactions = {}
    for index, reactant in enumerate(mapped_reaction.GetReactants()):
        unique_reaction = AllChem.ChemicalReaction()

        unique_reaction.AddReactantTemplate(reactant)
        for product in mapped_reaction.GetProducts():
            unique_reaction.AddProductTemplate(product)

        all_unique_reactions[f"reactantIdx_{index+1}"] = {
            "ID": reactant.GetProp("ID"),
            "SingleReactantReaction": _sort_reaction(unique_reaction),
        }

    return all_unique_reactions


def generate_rules(single_reactant_reaction: AllChem.ChemicalReaction):
    """
    Generates a dictionary of rules for a single reactant reaction based on the reacting atoms and the rings.

    Parameters
    ----------
    single_reactant_reaction: AllChem.ChemicalReaction
        The input chemical reaction with one reactant and one or more products.

    Returns
    -------
    reaction_rules: dict
        A dictionary of rules for the single reactant reaction, keyed by the reactant name and containing the reactant and product SMILES, the product name, and a sub-dictionary of rules keyed by the number of atoms to keep.

    Examples
    --------
    >>> reaction_smiles = "[H][C:1]([H])=[O:2]>>[O:2]=[C:1]([H])"
    >>> reaction = Chem.ReactionFromSmarts(reaction_smiles)
    >>> rules = generate_rules(reaction)
    >>> print(rules["reactant_1"]["ReactantMap"])
    >>> print(rules["reactant_1"]["ProductName"])
    >>> print(rules["reactant_1"]["ProductMap"])
    >>> print(rules["reactant_1"]["SingleReactantRules"][3])  # Rules for reactions with 3 atoms to keep
    """
    single_reactant_reaction.Initialize()

    reactant = single_reactant_reaction.GetReactantTemplate(0)
    reactant_name = reactant.GetProp("ID")

    main_product = single_reactant_reaction.GetProductTemplate(0)
    main_product_name = main_product.GetProp("ID")

    reacting_atoms = list(single_reactant_reaction.GetReactingAtoms()[0])

    matching_rings = []
    Chem.GetSSSR(reactant)
    rings = reactant.GetRingInfo()

    ## Initialize dictionaries to store data
    rules = {}
    reaction_rules = {}

    def _trim_reaction(single_reactant_reaction: AllChem.ChemicalReaction, atom_map_to_remove):
        """
        Trims a chemical reaction by removing specific atoms based on atom mapping.

        Parameters
        ----------
        single_reactant_reaction: AllChem.ChemicalReaction
            The input chemical reaction with one reactant and one or more products.
        atom_map_to_remove: list
            A list of atom map numbers corresponding to the atoms to be removed from the reaction.

        Returns
        -------
        trimmed_reaction: AllChem.ChemicalReaction
            The trimmed chemical reaction after removing the specified atoms.
        """
        single_reactant_reaction.Initialize()
        trimmed_reaction = AllChem.ChemicalReaction()

        def __adjust_mol_query(mol: Chem.Mol):
            params = Chem.AdjustQueryParameters()
            params.adjustDegree = False
            params.aromatizeIfPossible = True
            params.adjustRingChain = True
            params.adjustRingChainFlags = Chem.ADJUST_IGNOREDUMMIES
            params.adjustRingCount = True
            params.adjustSingleBondsBetweenAromaticAtoms = True
            params.adjustSingleBondsToDegreeOneNeighbors = True
            params.adjustConjugatedFiveRings = True

            fixed_mol = Chem.AdjustQueryProperties(
                mol,
                params,
            )

            for atom in fixed_mol.GetAtoms():
                if atom.GetIsAromatic():
                    atom.ExpandQuery(rdqueries.IsAromaticQueryAtom(False))

            return fixed_mol

        def __trim_mol(mol: Chem.Mol, atom_map_to_remove):
            fixed_mol = __adjust_mol_query(mol)
            editable_mol = Chem.EditableMol(fixed_mol)
            atoms_to_remove = [
                atom.GetIdx()
                for atom in mol.GetAtoms()
                if atom.GetAtomMapNum() in atom_map_to_remove
            ]

            for atom_index in sorted(atoms_to_remove, reverse=True):
                editable_mol.RemoveAtom(atom_index)

            return editable_mol.GetMol()

        for reactant in single_reactant_reaction.GetReactants():
            trimmed_reaction.AddReactantTemplate(
                __trim_mol(reactant, atom_map_to_remove)
            )
        for product in single_reactant_reaction.GetProducts():
            trimmed_reaction.AddProductTemplate(__trim_mol(product, atom_map_to_remove))
        return AllChem.ReactionToSmarts(trimmed_reaction)

    if len(reacting_atoms) < 1:
        rules = {}
    else:
        while (
            len(set(reacting_atoms)) <= 30
            and len(set(reacting_atoms)) < reactant.GetNumAtoms()
        ):
            for atom_index in set(reacting_atoms):
                atom = reactant.GetAtomWithIdx(atom_index)
                reacting_atoms.extend(
                    [
                        atom
                        for ring in rings.AtomRings()
                        for atom in ring
                        if atom_index in ring
                    ]
                )
                reacting_atoms.extend(
                    [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
                )
                atoms_to_keep = set(reacting_atoms)
                atom_map_to_remove = [
                    atom.GetAtomMapNum()
                    for atom in reactant.GetAtoms()
                    if atom.GetIdx() not in atoms_to_keep
                ]

            rules[len(atoms_to_keep)] = _trim_reaction(
                single_reactant_reaction, atom_map_to_remove
            )

    reaction_rules[reactant_name] = {
        "ReactantMap": dm.to_smiles(reactant),
        "ProductName": main_product_name,
        "ProductMap": dm.to_smiles(main_product),
        "SingleReactantRules": rules,
    }

    return reaction_rules


class Reaction(object):
    """
    A class to represent a chemical reaction with various attributes and methods.

    Attributes
    ----------
    SanitizedReaction: AllChem.ChemicalReaction
        The sanitized version of the input reaction, obtained by calling the SanitizeReaction function.
    MappedReaction: AllChem.ChemicalReaction
        The mapped version of the sanitized reaction, obtained by calling the MapReaction and SetReactionIds functions.
    ReversedReaction: AllChem.ChemicalReaction or None
        The reversed version of the mapped reaction, obtained by calling the ReverseReaction function, or None if the reversible argument is False.

    Methods
    -------
    __init__(reaction_smiles, reaction_ids, reversible, mapper)
        Initializes a REACTION object with the given arguments.

        Parameters
        -----------
            reaction_smiles:str
                The input chemical reaction in SMILES format.
            reaction_ids:str
                The IDs of the reactants and products of the input reaction, in the format of 'R1.R2.>>P1.P2.', where R1, R2, P1, P2 are the IDs.
            reversible:bool
                A flag to indicate whether to generate a reversed reaction or not. Default is False.
            mapper:str
                The name of the mapper to use for atom mapping. Either 'RXNMapper' or 'ReactionDecoder'. Default is 'ReactionDecoder'.

    Example
    -------
    >>> reaction_smiles = 'CC(=O)O.CCOC(=O)C>>CCOC(=O)CC.O'
    >>> reaction_ids = 'R1.R2>>P1.P2'
    >>> reaction = REACTION(reaction_smiles, reaction_ids, reversible=True, mapper='RXNMapper')
    >>> print(reaction.SanitizedReaction)
    [CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][O:7][C:8](=[O:9])[CH3:10]>><[CH3:1][CH2:6][O:7][C:8](=[O:9])[CH2:11][CH3:10].[OH:4][C:2](=[O:3])[H]
    >>> print(reaction.MappedReaction)
    [CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][O:7][C:8](=[O:9])[CH3:10]>><[CH3:1][CH2:6][O:7][C:8](=[O:9])[CH2:11][CH3:10].[OH:4][C:2](=[O:3])[H]
    >>> print(reaction.ReversedReaction)
    [CH3:1][CH2:6][O:7][C:8](=[O:9])[CH2:11][CH3:10].[OH:4][C:2](=[O:3])[H]>><[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][O:7][C:8](=[O:9])[CH3:10]
    """

    def __init__(
        self,
        reaction_smiles: str,
        reaction_ids: str,
        reversible: bool = False,
        mapper: str = "ReactionDecoder",
    ):
        self.SanitizedReaction = sanitize_reaction(
            AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
        )
        self.__mapped_reaction_raw = map_reaction(
            AllChem.ReactionToSmiles(self.SanitizedReaction), mapper=mapper
        )
        self.__mapped_reaction_sanitized = sanitize_reaction(self.__mapped_reaction_raw)
        self.MappedReaction = set_reaction_ids(
            self.SanitizedReaction, self.__mapped_reaction_sanitized, reaction_ids
        )
        if reversible == True:
            self.ReversedReaction = reverse_reaction(self.MappedReaction)
        if reversible == False:
            self.ReversedReaction = None