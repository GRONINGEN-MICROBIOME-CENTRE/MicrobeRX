import copy, subprocess
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdqueries, DataStructs,MACCSkeys

import datamol as dm

from rxnmapper import RXNMapper
rxn_mapper = RXNMapper()

import pandas as pd
import re

''' 
###################  REACTION TOOLKIT ###################
'''

from importlib_resources import files

REACTION_DECODER = files("microberx.RuleGenerator.bin").joinpath("RTD.jar")

''' 
###################  REACTION TOOLKIT ###################
'''

def SanitizeReaction (target_reaction:AllChem.ChemicalReaction) -> AllChem.ChemicalReaction:
    '''
    Sanitizes a chemical reaction by removing stereochemistry, replacing dummy atoms with carbon, and standardizing the molecules.

    Parameters:
        target_reaction (AllChem.ChemicalReaction): The input chemical reaction to be sanitized.

    Returns:
        fixed_reaction (AllChem.ChemicalReaction): The output chemical reaction after sanitization.
    '''
    def _replace_dummy_atoms(mol: Chem.Mol): ## Performs a molecular sanitization, removes stereochemistry and replaces dummy atoms
        dummy_query = Chem.MolFromSmiles('*')
        fixed_mol=AllChem.ReplaceSubstructs(mol,dummy_query,Chem.MolFromSmiles('[#6]'),replaceAll=True,useChirality=False)[0]
        fixed_mol = dm.fix_mol(fixed_mol)
        fixed_mol = dm.sanitize_mol(fixed_mol)
        fixed_mol = dm.standardize_mol(fixed_mol)
        Chem.RemoveStereochemistry(fixed_mol)
        AllChem.Compute2DCoords(fixed_mol)
        return fixed_mol


    fixed_reaction=AllChem.ChemicalReaction()
    
    for index, reactant in enumerate(target_reaction.GetReactants()):
        fixed_reactant=_replace_dummy_atoms(reactant)
        if fixed_reactant != None:
            fixed_reaction.AddReactantTemplate(fixed_reactant)
    for index, product in enumerate(target_reaction.GetProducts()):
        fixed_product=_replace_dummy_atoms(product)
        if fixed_product != None:
            fixed_reaction.AddProductTemplate(fixed_product)

    return fixed_reaction


def MapReaction (reaction_smiles:str=None, mapper:str='RXNMapper') -> AllChem.ChemicalReaction:
    '''
    Maps the atoms of a chemical reaction using either RXNMapper or ReactionDecoder.

    Parameters:
        reaction_smiles (str): The input chemical reaction in SMILES format.
        mapper (str): The name of the mapper to use. Either 'RXNMapper' or 'ReactionDecoder'. Default is 'RXNMapper'.

    Returns:
        mapped_reaction (AllChem.ChemicalReaction): The output rdkit chemical reaction with atom mapping.
    '''
    def _reaction_decoder(reaction_smiles:str):
        out=subprocess.call(['java', '-jar', REACTION_DECODER,'-Q', 'SMI', '-q', reaction_smiles, '-c', '-j', 'AAM', '-f', 'TEXT'],timeout=360)
        mapped_rxn=AllChem.ReactionFromRxnFile('ECBLAST_smiles_AAM.rxn')
        return mapped_rxn

    if mapper == 'RXNMapper':
        mapper=rxn_mapper.get_attention_guided_atom_maps([reaction_smiles])[0]

        mapped_reaction = AllChem.ReactionFromSmarts(mapper['mapped_rxn'],useSmiles=True)

    if mapper == 'ReactionDecoder':
        mapped_reaction= _reaction_decoder(reaction_smiles)

    return mapped_reaction

def SetReactionIds(reference_reaction:AllChem.ChemicalReaction, target_reaction:AllChem.ChemicalReaction, reaction_ids:str) -> AllChem.ChemicalReaction:
    '''
    Sets the IDs of the reactants and products of a target reaction based on their molecular formulas and a reference reaction.

    Parameters:
        reference_reaction (AllChem.ChemicalReaction): The reference chemical reaction that has the same reactants and products as the target reaction, but in a different order or orientation.
        target_reaction (AllChem.ChemicalReaction): The target chemical reaction that needs to have its IDs set.
        reaction_ids (str): The IDs of the reactants and products of the reference reaction, in the format of 'R1.R2.>>P1.P2.', where R1, R2, P1, P2 are the IDs.

    Returns:
        target_reaction (AllChem.ChemicalReaction): The target chemical reaction with its IDs set according to the reference reaction and the reaction_ids.
    '''
    reactants_ids=reaction_ids.split('>>')[0].split('.')
    products_ids=reaction_ids.split('>>')[1].split('.')
    
    reactants_formulas=[AllChem.CalcMolFormula(r) for r in reference_reaction.GetReactants()]
    products_formulas=[AllChem.CalcMolFormula(r) for r in reference_reaction.GetProducts()]

    reactants_map=dict(zip(reactants_formulas,reactants_ids))
    products_map=dict(zip(products_formulas,products_ids))
    
    for r in target_reaction.GetReactants():
        r.SetProp('ID',reactants_map[AllChem.CalcMolFormula(r)])

    for p in target_reaction.GetProducts():
        p.SetProp('ID',products_map[AllChem.CalcMolFormula(p)])
    
    return target_reaction
    

def ReverseReaction(reaction:AllChem.ChemicalReaction) -> AllChem.ChemicalReaction:
    '''
    ReverseReaction(reaction)
    Reverses a chemical reaction by swapping the reactants and products.

    Parameters:
        reaction (AllChem.ChemicalReaction): The input chemical reaction to be reversed.

    Returns:
        reversed_reaction (AllChem.ChemicalReaction): The output chemical reaction that is the reverse of the input reaction.
    '''
    reversed_reaction = AllChem.ChemicalReaction()
    
    for i, p in enumerate(reaction.GetProducts()):
        reversed_reaction.AddReactantTemplate(p)
    
    for i, r in enumerate(reaction.GetReactants()):
        reversed_reaction.AddProductTemplate(r)
    
    return reversed_reaction

''' 
###################  RULES TOOLKIT ###################
'''

def GenerateSingleReactantReactions(mapped_reaction:AllChem.ChemicalReaction) -> dict:
    '''
    Generates a dictionary of single reactant reactions from a mapped reaction.

    Parameters:
        mapped_reaction (AllChem.ChemicalReaction): The input chemical reaction with atom mapping.

    Returns:
        all_unique_reactions (dict): A dictionary of single reactant reactions, keyed by the reactant index and containing the reactant ID and the single reactant reaction object.

    Helper function:
        _sort_reaction(single_reactant_reaction)
            Sorts the products of a single reactant reaction based on the number of mapped atoms they share with the reactant.

            Parameters:
                single_reactant_reaction (AllChem.ChemicalReaction): The input chemical reaction with one reactant and one or more products.

            Returns:
                sorted_reaction (AllChem.ChemicalReaction): The output chemical reaction with the same reactant and the products sorted in descending order of mapped atom overlap.
    '''
        
    def _sort_reaction(single_reactant_reaction:AllChem.ChemicalReaction) -> AllChem.ChemicalReaction:

        reactant_atom_map=[atom.GetAtomMapNum() for atom in single_reactant_reaction.GetReactantTemplate(0).GetAtoms()]

        if len(reactant_atom_map)>0:
            sorted_reaction=AllChem.ChemicalReaction()
            [sorted_reaction.AddReactantTemplate(reactant) for reactant in single_reactant_reaction.GetReactants()]
            products={}
            for product in single_reactant_reaction.GetProducts():
                product_atom_map=[atom.GetAtomMapNum() for atom in product.GetAtoms() if atom.GetAtomMapNum()!=0]
                products[product]=len(set(product_atom_map) & set(reactant_atom_map))

            sorted_products=dict(sorted(products.items(), key = lambda x: x[1], reverse = True))

            [sorted_reaction.AddProductTemplate(product) for product in sorted_products.keys()]

        return sorted_reaction


    all_unique_reactions={}
    for index,reactant in enumerate(mapped_reaction.GetReactants()):
        unique_reaction=AllChem.ChemicalReaction()

        unique_reaction.AddReactantTemplate(reactant)
        for product in mapped_reaction.GetProducts():
            unique_reaction.AddProductTemplate(product)

        all_unique_reactions[f"reactantIdx_{index+1}"]={'ID':reactant.GetProp('ID'),'SingleReactantReaction':_sort_reaction(unique_reaction)}

    return all_unique_reactions

def GenerateRules(single_reactant_reaction:AllChem.ChemicalReaction) -> dict:
    '''
    GenerateRules(single_reactant_reaction)
    Generates a dictionary of rules for a single reactant reaction based on the reacting atoms and the rings.

    Parameters:
        single_reactant_reaction (AllChem.ChemicalReaction): The input chemical reaction with one reactant and one or more products.

    Returns:
        reaction_rules (dict): A dictionary of rules for the single reactant reaction, keyed by the reactant name and containing the reactant and product SMILES, the product name, and a sub-dictionary of rules keyed by the number of atoms to keep.

    Helper function:
        _trim_reaction(single_reactant_reaction, atom_map_to_remove)
            Trims a single reactant reaction by removing the atoms with the specified atom map numbers.

            Parameters:
                single_reactant_reaction (AllChem.ChemicalReaction): The input chemical reaction with one reactant and one or more products.
                atom_map_to_remove (list): A list of atom map numbers to remove from the reaction.

            Returns:
                trimmed_reaction (str): The output chemical reaction in SMARTS format after trimming.

            Helper function:
                __adjust_mol_query(mol)
                    Adjusts the query properties of a molecule to make it more generic.

                    Parameters:
                        mol (Chem.Mol): The input molecule to be adjusted.

                    Returns:
                        fixed_mol (Chem.Mol): The output molecule after adjusting the query properties.

                __trim_mol(mol, atom_map_to_remove)
                    Trims a molecule by removing the atoms with the specified atom map numbers.

                    Parameters:
                        mol (Chem.Mol): The input molecule to be trimmed.
                        atom_map_to_remove (list): A list of atom map numbers to remove from the molecule.

                    Returns:
                        trimmed_mol (Chem.Mol): The output molecule after trimming.
    '''
    single_reactant_reaction.Initialize()
    
    reactant = single_reactant_reaction.GetReactantTemplate(0)
    reactant_name = reactant.GetProp('ID')
    
    main_product = single_reactant_reaction.GetProductTemplate(0)
    main_product_name = main_product.GetProp('ID')
    
    reacting_atoms=list(single_reactant_reaction.GetReactingAtoms()[0])   

    matching_rings=[]
    Chem.GetSSSR(reactant)
    rings=reactant.GetRingInfo()
    
    ## Initialize dictionaries to store data
    rules = {}
    reaction_rules = {}

    def _trim_reaction(single_reactant_reaction: AllChem.ChemicalReaction, atom_map_to_remove) -> AllChem.ChemicalReaction:            
        single_reactant_reaction.Initialize()
        trimmed_reaction= AllChem.ChemicalReaction()

        def __adjust_mol_query(mol:Chem.Mol):

            params = Chem.AdjustQueryParameters()
            params.adjustDegree=False
            params.aromatizeIfPossible=True
            params.adjustRingChain=True
            params.adjustRingChainFlags = Chem.ADJUST_IGNOREDUMMIES
            params.adjustRingCount=True
            params.adjustSingleBondsBetweenAromaticAtoms=True
            params.adjustSingleBondsToDegreeOneNeighbors=True
            params.adjustConjugatedFiveRings=True

            fixed_mol = Chem.AdjustQueryProperties(mol,params,)

            for atom in fixed_mol.GetAtoms():
                if atom.GetIsAromatic(): 
                    atom.ExpandQuery(rdqueries.IsAromaticQueryAtom(False))

            return fixed_mol

        def __trim_mol(mol:Chem.Mol,atom_map_to_remove):

            fixed_mol=__adjust_mol_query(mol)
            editable_mol = Chem.EditableMol(fixed_mol)
            atoms_to_remove= [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() in atom_map_to_remove]

            for atom_index in sorted(atoms_to_remove, reverse=True):
                editable_mol.RemoveAtom(atom_index)

            return editable_mol.GetMol()


        for reactant in single_reactant_reaction.GetReactants():
            trimmed_reaction.AddReactantTemplate(__trim_mol(reactant,atom_map_to_remove))
        for product in single_reactant_reaction.GetProducts():
            trimmed_reaction.AddProductTemplate(__trim_mol(product,atom_map_to_remove))
        return AllChem.ReactionToSmarts(trimmed_reaction)
    
    if len(reacting_atoms)<1:
        rules={}
    else:
    
        while len(set(reacting_atoms))<=30 and len(set(reacting_atoms))<reactant.GetNumAtoms():
                for atom_index in set(reacting_atoms):
                    atom=reactant.GetAtomWithIdx(atom_index)
                    reacting_atoms.extend([atom for ring in rings.AtomRings() for atom in ring if atom_index in ring])
                    reacting_atoms.extend([neighbor.GetIdx() for neighbor in atom.GetNeighbors()])
                    atoms_to_keep=(set(reacting_atoms))
                    atom_map_to_remove=[atom.GetAtomMapNum() for atom in reactant.GetAtoms() if atom.GetIdx() not in atoms_to_keep] 

                rules[len(atoms_to_keep)]=_trim_reaction(single_reactant_reaction,atom_map_to_remove)
    
    reaction_rules[reactant_name]={'ReactantMap':dm.to_smiles(reactant), 'ProductName': main_product_name,'ProductMap':dm.to_smiles(main_product),'SingleReactantRules':rules}

    return reaction_rules

'''
###################  CLASSES  ###################
'''

class REACTION(object):
    '''
    A class to represent a chemical reaction with various attributes and methods.

    Attributes:
        SanitizedReaction (AllChem.ChemicalReaction): The sanitized version of the input reaction, obtained by calling the SanitizeReaction function.
        MappedReaction (AllChem.ChemicalReaction): The mapped version of the sanitized reaction, obtained by calling the MapReaction and SetReactionIds functions.
        ReversedReaction (AllChem.ChemicalReaction or None): The reversed version of the mapped reaction, obtained by calling the ReverseReaction function, or None if the reversible argument is False.

    Methods:
        __init__(reaction_smiles, reaction_ids, reversible, mapper)
            Initializes a REACTION object with the given arguments.

            Parameters:
                reaction_smiles (str): The input chemical reaction in SMILES format.
                reaction_ids (str): The IDs of the reactants and products of the input reaction, in the format of 'R1.R2.>>P1.P2.', where R1, R2, P1, P2 are the IDs.
                reversible (bool): A flag to indicate whether to generate a reversed reaction or not. Default is False.
                mapper (str): The name of the mapper to use for atom mapping. Either 'RXNMapper' or 'ReactionDecoder'. Default is 'ReactionDecoder'.

    Example:
        >>> from rdkit import Chem
        >>> from rdkit.Chem import AllChem
        >>> import dm # a module for molecular sanitization and standardization
        >>> import rxn_mapper # a module for attention-guided atom mapping
        >>> import subprocess # a module for calling external commands
        >>> reaction_smiles = 'CC(=O)O.CCOC(=O)C>>CCOC(=O)CC.O'
        >>> reaction_ids = 'R1.R2>>P1.P2'
        >>> reaction = REACTION(reaction_smiles, reaction_ids, reversible=True, mapper='RXNMapper')
        >>> print(reaction.SanitizedReaction)
        [CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][O:7][C:8](=[O:9])[CH3:10]>><[CH3:1][CH2:6][O:7][C:8](=[O:9])[CH2:11][CH3:10].[OH:4][C:2](=[O:3])[H]
        >>> print(reaction.MappedReaction)
        [CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][O:7][C:8](=[O:9])[CH3:10]>><[CH3:1][CH2:6][O:7][C:8](=[O:9])[CH2:11][CH3:10].[OH:4][C:2](=[O:3])[H]
        >>> print(reaction.ReversedReaction)
        [CH3:1][CH2:6][O:7][C:8](=[O:9])[CH2:11][CH3:10].[OH:4][C:2](=[O:3])[H]>><[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][CH2:6][O:7][C:8](=[O:9])[CH3:10]
        '''
    def __init__(self, reaction_smiles:str, reaction_ids:str, reversible:bool = False, mapper:str = 'ReactionDecoder'):
        self.SanitizedReaction = SanitizeReaction(AllChem.ReactionFromSmarts(reaction_smiles,useSmiles=True))
        self.__mapped_reaction_raw = MapReaction(AllChem.ReactionToSmiles(self.SanitizedReaction),mapper=mapper)
        self.__mapped_reaction_sanitized = SanitizeReaction(self.__mapped_reaction_raw)
        self.MappedReaction = SetReactionIds(self.SanitizedReaction,self.__mapped_reaction_sanitized,reaction_ids)
        if reversible==True:
            self.ReversedReaction = ReverseReaction(self.MappedReaction)
        if reversible==False:
            self.ReversedReaction = None

class PARSER(object):
    '''
    A class to parse a chemical reaction from a string format into a dictionary format.

    Attributes:
        reaction (str): The input chemical reaction in a string format, using names or symbols of compounds and stoichiometry coefficients.
        compounds_map (dict): A dictionary that maps the names or symbols of compounds to their SMILES strings.
        reaction_dict (dict): The output chemical reaction in a dictionary format, containing the reversibility, the reactants and products with their stoichiometry and SMILES, and the reaction SMILES and names.

    Methods:
        __init__(reaction, compounds_map)
            Initializes a PARSER object with the given arguments.

            Parameters:
                reaction (str): The input chemical reaction in a string format, using names or symbols of compounds and stoichiometry coefficients.
                compounds_map (dict): A dictionary that maps the names or symbols of compounds to their SMILES strings.

        decompose_reaction()
            Decomposes the input reaction into its components and stores them in the reaction_dict attribute.

            No parameters or returns.

    Example:
        >>> import re # a module for regular expressions
        >>> reaction = '2 H2 + O2 --> 2 H2O'
        >>> compounds_map = {'H2': '[HH]', 'O2': 'O=O', 'H2O': 'O'}
        >>> parser = PARSER(reaction, compounds_map)
        >>> parser.decompose_reaction()
        >>> print(parser.reaction_dict)
        {'Reversible': False, 'LEFT': {'H2': {'stoichiometry': 2.0, 'smiles': '[HH]'}, 'O2': {'stoichiometry': 1.0, 'smiles': 'O=O'}}, 'RIGHT': {'H2O': {'stoichiometry': 2.0, 'smiles': 'O'}}, 'ReactionSmiles': '[HH].[HH].O=O>>O.O', 'ReactionNames': 'H2.H2.O2>>H2O.H2O'}
    '''
    def __init__(self,reaction:str=None,compounds_map:dict=None):
        self.reaction=reaction
        self.compounds_map=compounds_map

    def decompose_reaction(self):
        if ' = ' in self.reaction:
            self.reaction = re.sub("@\w+", "", self.reaction)
            reactants, products = self.reaction.split(" = ")
            Reversible=True
        
        if '<=>' in self.reaction:
            self.reaction = re.sub(r"\[[a-z]\]", "", self.reaction)
            reactants, products = self.reaction.split('<=>')
            Reversible=True
        
        if '-->' in self.reaction:
            self.reaction = re.sub(r"\[[a-z]\]", "", self.reaction)
            reactants, products = self.reaction.split('-->')
            Reversible=False
        
        
        # Split the reactants and products into individual compounds
        reactants = reactants.split('+')
        products = products.split('+')

        # Create empty dictionaries to store the reactants and products with their stoichiometry
        Left_side = {}
        Right_side = {}
        smiles_scheme=""
        names_scheme=""
        
        for r in reactants:
            r = r.strip()
            if len(r.split())>1:
                stoich, name = r.split()
            else:
                name=r
                stoich=1
            
            smiles=self.compounds_map[name]
            Left_side[name] = {'stoichiometry':float(stoich),'smiles':smiles}
            smiles_scheme += (smiles + ".")* int(float(stoich))
            names_scheme += (name + ".")* int(float(stoich))
        
        smiles_scheme = smiles_scheme[:-1]
        names_scheme = names_scheme[:-1]
        
        smiles_scheme += ">>"
        names_scheme += ">>"
        for p in products:
            p = p.strip()
            if len(p.split())>1:
                stoich, name = p.split()
            else:
                name=p
                stoich=1
            
            smiles=self.compounds_map[name]
            Right_side[name] = {'stoichiometry':float(stoich),'smiles':smiles}
            smiles_scheme += (smiles + ".")* int(float(stoich))
            names_scheme += (name + ".")* int(float(stoich))

        self.reaction_dict={'Reversible':Reversible,'LEFT':Left_side,'RIGHT':Right_side,'ReactionSmiles':smiles_scheme[:-1],'ReactionNames':names_scheme[:-1]}