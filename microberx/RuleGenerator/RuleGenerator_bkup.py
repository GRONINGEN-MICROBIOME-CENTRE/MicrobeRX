import copy
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdqueries, DataStructs,MACCSkeys

from rxnmapper import RXNMapper
rxn_mapper = RXNMapper()

import pandas as pd


class Reaction(object):
    def __init__(self,reaction_smiles:str, reaction_ids:str, id_key:str='ID'):
        self.__reaction_with_dummies = AllChem.ReactionFromSmarts(reaction_smiles,useSmiles=True)
        self.reaction_ids = reaction_ids
        self.reactant_ids=self.reaction_ids.split('>>')[0].split('.')
        self.product_ids=self.reaction_ids.split('>>')[1].split('.')
        self.id_key=id_key
        self.fixed_reaction = self.__set_ids_to_reaction_and_adjust_agents()
        self.fixed_reaction_smiles = AllChem.ReactionToSmiles(self.fixed_reaction)
        
        self.mapped_reaction = self.__map_reaction()
        self.mapped_reaction_smarts = AllChem.ReactionToSmarts (self.mapped_reaction)
        
    def __set_ids_to_reaction_and_adjust_agents (self):
        
        def ___replace_dummy_atoms(mol: Chem.Mol): ## Performs a molecular sanitization, removes stereochemistry and replaces dummy atoms
            Chem.SanitizeMol(mol)
            Chem.RemoveStereochemistry(mol)
            dummy_query = Chem.MolFromSmiles('*')
            fixed_mol=AllChem.ReplaceSubstructs(mol,dummy_query,Chem.MolFromSmiles('[#6]'),replaceAll=True,useChirality=False)[0]
            #fixed_mol= ___neutralize_atoms(fixed_mol)
            AllChem.Compute2DCoords(fixed_mol)
            return fixed_mol
        
        def ___neutralize_atoms(mol:Chem.Mol):
            charged_atom_smarts = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
            atom_matches = mol.GetSubstructMatches(charged_atom_smarts)
            atom_matches_list = [atom_list[0] for atom_list in atom_matches]
            if len(atom_matches_list) > 0:
                for atom_index in atom_matches_list:
                    atom = mol.GetAtomWithIdx(atom_index)
                    atom_charge = atom.GetFormalCharge()
                    atom_hcount = atom.GetTotalNumHs()
                    atom.SetFormalCharge(0)
                    atom.SetNumExplicitHs(atom_hcount - atom_charge)
                    atom.UpdatePropertyCache()
        
        fixed_reaction=AllChem.ChemicalReaction()
        for index, reactant in enumerate(self.__reaction_with_dummies.GetReactants()):
            fixed_reactant=___replace_dummy_atoms(reactant)
            if fixed_reactant != None:
                fixed_reactant.SetProp(self.id_key,self.reactant_ids[index])
                fixed_reaction.AddReactantTemplate(fixed_reactant)
        for index, product in enumerate(self.__reaction_with_dummies.GetProducts()):
            fixed_product=___replace_dummy_atoms(product)
            if fixed_product != None:
                fixed_product.SetProp(self.id_key,self.product_ids[index])
                fixed_reaction.AddProductTemplate(fixed_product)
        
        return fixed_reaction
    
    def __map_reaction (self):
        
        def ___reset_molecules_ids(reference_mols:list, target_mols: list):
            [t.SetProp(self.id_key,'No_ID') for t in target_mols]
            for reference in reference_mols:
                Chem.GetSSSR(reference)
                reference_fps = MACCSkeys.GenMACCSKeys(reference)
                for target in target_mols:
                    Chem.GetSSSR(target)
                    target_fps = MACCSkeys.GenMACCSKeys(target)
                    similarity=DataStructs.FingerprintSimilarity(reference_fps,target_fps)
                    if similarity == 1:
                        target.SetProp(self.id_key, reference.GetProp(self.id_key))
                        break
                        
        mapper=rxn_mapper.get_attention_guided_atom_maps([self.fixed_reaction_smiles])[0]

        mapped_reaction = AllChem.ReactionFromSmarts(mapper['mapped_rxn'],useSmiles=True)
        
        ___reset_molecules_ids (self.fixed_reaction.GetReactants(), mapped_reaction.GetReactants())
        ___reset_molecules_ids (self.fixed_reaction.GetProducts(), mapped_reaction.GetProducts())
        
        return mapped_reaction
        
            
    def __rename_reaction_agents_after_mapping (self):
        [tarMol.SetProp(self.id_key, 'X') for tarMol in tarMols]
        for refMol in refMols:
            refFPS = Chem.RDKFingerprint(refMol, useHs=True)
            for tarMol in tarMols:
                tarFPS = Chem.RDKFingerprint(tarMol, useHs=True)
                sim = DataStructs.FingerprintSimilarity(refFPS, tarFPS)
                if sim == 1:
                    tarMol.SetProp(self.id_key, refMol.GetProp(self.id_key))
                    break
                    

class Rules(object):
    def __init__(self, mapped_reaction:AllChem.ChemicalReaction,id_key:str='ID'):
        self.mapped_reaction=mapped_reaction
        self.id_key=id_key
        self.single_reactant_reactions = self.__generate_single_reactant_reactions()
        self.single_reactant_rules={}
        for key in self.single_reactant_reactions:
            rules=self.__generate_rules(self.single_reactant_reactions[key])
            self.single_reactant_rules[key]=rules
    
    def __generate_single_reactant_reactions(self):
        
        def ___sort_reaction(single_reactant_reaction:AllChem.ChemicalReaction) -> AllChem.ChemicalReaction:
            
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
        for reactant in self.mapped_reaction.GetReactants():
            unique_reaction=AllChem.ChemicalReaction()
            
            unique_reaction.AddReactantTemplate(reactant)
            for product in self.mapped_reaction.GetProducts():
                unique_reaction.AddProductTemplate(product)
        
            all_unique_reactions[reactant.GetProp(self.id_key)]=___sort_reaction(unique_reaction)
        
        return all_unique_reactions
    
    def __generate_rules(self,single_reactant_reaction:AllChem.ChemicalReaction):
        rules={}
        single_reactant_reaction.Initialize()
        reactant=single_reactant_reaction.GetReactantTemplate(0)
        reacting_atoms=list(single_reactant_reaction.GetReactingAtoms()[0])   
        
        matching_rings=[]
        Chem.GetSSSR(reactant)
        rings=reactant.GetRingInfo()
        
        def ___trim_reaction(single_reactant_reaction: AllChem.ChemicalReaction, atom_map_to_remove) -> AllChem.ChemicalReaction:
            single_reactant_reaction.Initialize()
            trimmed_reaction= AllChem.ChemicalReaction()
            
            def ____adjust_mol_query(mol:Chem.Mol):

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
            
            def ____trim_mol(mol:Chem.Mol,atom_map_to_remove):
                
                fixed_mol=____adjust_mol_query(mol)
                editable_mol = Chem.EditableMol(fixed_mol)
                atoms_to_remove= [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() in atom_map_to_remove]
                
                for atom_index in sorted(atoms_to_remove, reverse=True):
                    editable_mol.RemoveAtom(atom_index)
                
                return editable_mol.GetMol()

            for reactant in single_reactant_reaction.GetReactants():
                trimmed_reaction.AddReactantTemplate(____trim_mol(reactant,atom_map_to_remove))
            for product in single_reactant_reaction.GetProducts():
                trimmed_reaction.AddProductTemplate(____trim_mol(product,atom_map_to_remove))

            return AllChem.ReactionToSmarts(trimmed_reaction)
        
        while len(set(reacting_atoms))<30 and len(set(reacting_atoms))<reactant.GetNumAtoms():
                for atom_index in set(reacting_atoms):
                    atom=reactant.GetAtomWithIdx(atom_index)
                    reacting_atoms.extend([atom for ring in rings.AtomRings() for atom in ring if atom_index in ring])
                    reacting_atoms.extend([neighbor.GetIdx() for neighbor in atom.GetNeighbors()])
                    atoms_to_keep=(set(reacting_atoms))
                    atom_map_to_remove=[atom.GetAtomMapNum() for atom in reactant.GetAtoms() if atom.GetIdx() not in atoms_to_keep]

                rules[len(atoms_to_keep)]=___trim_reaction(single_reactant_reaction,atom_map_to_remove)

        return rules