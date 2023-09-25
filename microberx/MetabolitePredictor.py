"""
This is a module that provides the main classes to predicted metabolites using MicrobeRX.

The module contains the following classes:

- MetabolitePredictor:  A class for predicting metabolites using reaction rules.

- RunPredictionRule: A class for predicting metabolites based on a single reaction rule.
"""

__all__ = ["MetabolitePredictor", "RunPredictionRule"]

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, DataStructs

import itertools, copy
import pandas as pd

import datamol as dm

from tqdm.notebook import tqdm

rdkit.RDLogger.DisableLog("rdApp.*")


class MetabolitePredictor:
    """
    A class for predicting metabolites using reaction rules.

    Parameters
    ----------
    rules_table : str
        The path to a table containing reaction rules and associated information.

    Attributes
    ----------
    predicted_metabolites : pd.DataFrame
        A DataFrame to store predicted metabolites and associated information.
    query : Chem.Mol
        The query molecule for metabolite prediction.
    query_name : str
        The name associated with the query molecule.
    query_atoms_num : int
        The number of heavy atoms in the query molecule.
    reacting_atoms : list
        A list to store reacting atom indices.
    reacting_atoms_in_unique_metabolites : list
        A list to store reacting atom indices in unique metabolites.

    Methods
    -------
    run_prediction(query: Chem.Mol, name: str = "metabolite")
        Run the metabolite prediction using the provided query molecule and name.
    """
    
    def __init__(self, rules_table: pd.DataFrame):
        """
        Initialize a MetabolitePredictor instance.

        Parameters
        ----------
        rules_table : pd.DataFrame
            The path to a table containing reaction rules and associated information.
        """
        self.rules_table = rules_table

    def run_prediction(self, query: Chem.Mol, name: str = "metabolite"):
        """
        This method performs metabolite prediction using reaction rules from the rules table. It calculates confidence scores for predicted metabolites and stores the results in the 'predicted_metabolites' attribute.

        Parameters
        ----------
        query : Chem.Mol
            The query molecule for metabolite prediction.
        name : str, optional
            The name associated with the query molecule. Default is "metabolite".

        Returns
        -------
        None
        
        Example
        -------
        >>> predictor = MetabolitePredictor(rules_table)
        >>> query_molecule = Chem.MolFromSmiles("CC(=O)O")
        >>> predictor.run_prediction(query_molecule, "acetate")
        >>> predicted_metabolites_df = predictor.predicted_metabolites
        """
        self.predicted_metabolites = pd.DataFrame()
        self.query = query
        self.query_name = name
        self.query_atoms_num = query.GetNumHeavyAtoms()
        self.reacting_atoms = []
        self.reacting_atoms_in_unique_metabolites = []

        for index in tqdm(self.rules_table.index):
            # try:
            Predictor = RunPredictionRule(
                query=self.query,
                rule_smarts=self.rules_table.rule[index],
                real_product=self.rules_table.product_map[index],
                real_substrate=self.rules_table.substrate_map[index],
            )
            Predictor.predict()
            results = Predictor.unique_products
            if results:
                tmp_table = pd.DataFrame(results.values())
                tmp_table["reaction_id"] = self.rules_table["reaction_id"][index]
                tmp_table["substrate"] = self.rules_table.substrate[index]
                tmp_table["product"] = self.rules_table["product"][index]
                tmp_table["num_atoms"] = self.rules_table.num_atoms[index]
                tmp_table["reacting_atoms_efficiency"] = round(
                    self.rules_table.num_atoms[index] / self.query_atoms_num, 3
                )
                self.predicted_metabolites = pd.concat(
                    [self.predicted_metabolites, tmp_table], ignore_index=True
                )
        # except Exception:
        # pass

        self.predicted_metabolites["confidence_score"] = round(
            self.predicted_metabolites.similarity_substrates
            + self.predicted_metabolites.similarity_products
            + self.predicted_metabolites.reacting_atoms_efficiency,
            3,
        )
        self.predicted_metabolites.sort_values(
            by="confidence_score", ascending=False, ignore_index=True, inplace=True
        )
        self.predicted_metabolites.drop_duplicates(
            subset=["reaction_id", "substrate", "product", "main_product_smiles"],
            keep="first",
            inplace=True,
            ignore_index=True,
        )
        # self.predicted_metabolites = self.predicted_metabolites.merge(self.evidences_database[['reaction_id','ec','reference_db','reaction_name','gene_name','uniprot','entrez']],on='reaction_id',how='left')

        un_mets = self.predicted_metabolites.drop_duplicates(
            subset=["main_product_smiles"], keep="first", ignore_index=True
        )
        self.unique_metabolites = copy.deepcopy(un_mets)
        metabolite_rank = {
            self.unique_metabolites.main_product_smiles[i]: f"{self.query_name}_{i+1}"
            for i in self.unique_metabolites.index
        }

        self.unique_metabolites[
            "metabolite_id"
        ] = self.unique_metabolites.main_product_smiles.map(metabolite_rank)
        self.predicted_metabolites[
            "metabolite_id"
        ] = self.predicted_metabolites.main_product_smiles.map(metabolite_rank)


class RunPredictionRule:
    """
    A class for predicting reactions based on reaction rules.

    Parameters
    ----------
    query : Chem.Mol
        The query molecule for prediction.
    rule_smarts : str
        The reaction rule in SMARTS format.
    real_product : str
        The SMILES representation of the real product.
    real_substrate : str
        The SMILES representation of the real substrate.

    Attributes
    ----------
    query : Chem.Mol
        The query molecule for prediction.
    rule_smarts : str
        The reaction rule in SMARTS format.
    reaction : AllChem.Reaction
        The reaction object created from the reaction rule.
    real_product : Chem.Mol
        The real product molecule.
    real_substrate : Chem.Mol
        The real substrate molecule.
    unique_products : dict
        A dictionary to store information about unique predicted products.

    Methods
    -------
    predict()
        Predict reaction products and populate the 'unique_products' attribute.

    """
    def __init__(self, query: Chem.Mol, rule_smarts: str, real_product: str, real_substrate: str):
        """
        Initialize a RunPredictionRule instance.

        Parameters
        ----------
        query : Chem.Mol
            The query molecule for prediction.
        rule_smarts : str
            The reaction rule in SMARTS format.
        real_product : str
            The SMILES representation of the real product.
        real_substrate : str
            The SMILES representation of the real substrate.
        """
        self.query = query
        self.rule_smarts = rule_smarts
        self.reaction = AllChem.ReactionFromSmarts(self.rule_smarts, useSmiles=False)
        self.reaction.Initialize()
        self.real_product = Chem.MolFromSmiles(real_product)
        self.real_substrate = Chem.MolFromSmiles(real_substrate)

    def predict(self):
        """
        This method predicts reaction products based on the provided query molecule, reaction rule, real product, and real substrate. It calculates various properties and stores the results in the 'unique_products' dictionary.

        Returns
        -------
        None
        """
        self._main_reactant = self.reaction.GetReactantTemplate(0)
        num_substructure_matches = len(
            self.query.GetSubstructMatches(self._main_reactant)
        )
        self.unique_products = {}
        if num_substructure_matches > 0:
            self._atom_map = [
                a.GetAtomMapNum()
                for a in self._main_reactant.GetAtoms()
                if a.GetAtomMapNum() != 0
            ]

            results = self.reaction.RunReactants([self.query], maxProducts=10)

            for products in results:
                try:
                    main_product = products[0]
                    fixed_mol = dm.fix_mol(main_product)
                    fixed_mol = dm.sanitize_mol(fixed_mol)
                    fixed_mol = dm.standardize_mol(fixed_mol)
                    smi = Chem.MolToSmiles(
                        fixed_mol, isomericSmiles=False, canonical=True
                    )

                    (
                        main_product_atom_indexes,
                        query_atoms_indexes,
                    ) = self.__get_product_atom_indexes(main_product)

                    real_product_atom_indexes = [
                        a.GetIdx()
                        for a in self.real_product.GetAtoms()
                        if a.GetAtomMapNum() in self._atom_map
                    ]
                    real_substrate_atom_indexes = [
                        a.GetIdx()
                        for a in self.real_substrate.GetAtoms()
                        if a.GetAtomMapNum() in self._atom_map
                    ]

                    query_substructure = self.__get_mol_substructure(
                        self.query, query_atoms_indexes
                    )
                    real_substrate_substructure = self.__get_mol_substructure(
                        self.real_substrate, real_substrate_atom_indexes
                    )
                    main_product_substructure = self.__get_mol_substructure(
                        main_product, main_product_atom_indexes
                    )
                    real_product_substructure = self.__get_mol_substructure(
                        self.real_product, real_product_atom_indexes
                    )

                    self.unique_products[smi] = {
                        "main_product_smiles": smi,
                        "secondary_products_smiles": ".".join(
                            [Chem.MolToSmiles(p, canonical=True) for p in products[1:]]
                        ),
                        "similarity_substrates": self.__compute_similarity(
                            self.query, self.real_substrate, fingerPrintModel="rdkit"
                        ),
                        #'similarity_substrates_substructure':self.__compute_similarity(query_substructure,real_substrate_substructure,fingerPrintModel='rdkit'),
                        "similarity_products": self.__compute_similarity(
                            main_product, self.real_product, fingerPrintModel="rdkit"
                        ),
                        #'similarity_products_substructure':self.__compute_similarity(main_product_substructure,real_product_substructure,fingerPrintModel='rdkit'),
                        "reacting_atoms_in_query": query_atoms_indexes,
                    }
                except Exception as e:
                    pass

    def __get_product_atom_indexes(self, product: Chem.Mol):
        product_atoms_indexes = []
        query_atoms_indexes = []
        for atom in product.GetAtoms():
            atom_properties = atom.GetPropsAsDict()
            if set(["old_mapno", "react_atom_idx"]).issubset(
                set(atom_properties.keys())
            ):
                product_atoms_indexes.append(atom.GetIdx())
                query_atoms_indexes.append(int(atom.GetProp("react_atom_idx")))
        return product_atoms_indexes, query_atoms_indexes

    def __get_mol_substructure(self, mol_target: Chem.Mol, atom_indexes_to_keep):
        """
        Get the reaction substructure from the input molecule.
        :return: The reaction substructure
        """
        atom_indexes_to_keep = self.__get_substructure_neighbors(
            mol_target, atom_indexes_to_keep
        )
        atom_indexes_to_remove = [
            atom.GetIdx()
            for atom in mol_target.GetAtoms()
            if atom.GetIdx() not in atom_indexes_to_keep
        ]
        substructure = self.__trim_mol(mol_target, atom_indexes_to_remove)

        return Chem.MolFromSmiles(Chem.MolToSmiles(substructure))

    def __get_substructure_neighbors(self, mol_target: Chem.Mol, atom_indexes_to_keep):
        """
        Get the substructure neighbors of the input molecule.
        :return: The substructure neighbors
        """

        neighbor_atoms = []
        rings = mol_target.GetRingInfo()
        for atom_index in atom_indexes_to_keep:
            atom = mol_target.GetAtomWithIdx(atom_index)
            neighbor_atoms.extend(
                [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
            )
            neighbor_atoms.extend(
                [
                    neighbor_neighbor.GetIdx()
                    for neighbor in atom.GetNeighbors()
                    for neighbor_neighbor in neighbor.GetNeighbors()
                ]
            )
            neighbor_atoms.extend(
                [a for ring in rings.AtomRings() for a in ring if atom_index in ring]
            )

        unique_atom_pairs = list(itertools.combinations(set(neighbor_atoms), r=2))
        atom_path = []
        for i, j in unique_atom_pairs:
            atom_path += Chem.GetShortestPath(mol_target, i, j)

        return list(set(atom_path))

    def __trim_mol(self, mol_target: Chem.Mol, atom_indexes_to_remove):
        """
        Trim the input molecule to get the desired substructure.
        :param atoms_to_remove: The atoms to remove from the input molecule
        :return: The trimmed molecule
        """
        mol_target = copy.deepcopy(mol_target)
        Chem.SanitizeMol(mol_target)
        Chem.Kekulize(mol_target, clearAromaticFlags=True)

        editable_mol = Chem.EditableMol(mol_target)
        for atom_id in sorted(atom_indexes_to_remove, reverse=True):
            editable_mol.RemoveAtom(atom_id)

        return editable_mol.GetMol()

    def __compute_similarity(
        self, mol1: Chem.Mol, mol2: Chem.Mol, fingerPrintModel: str = "maccskey"
    ):
        if fingerPrintModel == "rdkit":
            similarity = DataStructs.FingerprintSimilarity(
                Chem.RDKFingerprint(mol1), Chem.RDKFingerprint(mol2)
            )
        if fingerPrintModel == "maccskey":
            similarity = DataStructs.FingerprintSimilarity(
                MACCSkeys.GenMACCSKeys(mol1), MACCSkeys.GenMACCSKeys(mol2)
            )

        return round(similarity, 3)
