import os
import multiprocessing
import argparse
import pandas as pd
from microberx import MetabolitePredictor
from rdkit import Chem
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")

"""
MicrobeRX Predictions Script
=============================

This script allows users to perform predictions using the MicrobeRX library either with a single smiles string or a file containing multiple smiles strings. The predictions can be customized with various parameters such as the biosystem, cutoff value, number of processing cores, and output directory.

Usage
-----
To run the script, use the following command format:

```bash
python microberx.py --smiles <smiles_string> --query_name <query_name> --biosystem <biosystem> --cutoff <cutoff_value> --out <output_directory>

python microberx.py --file <input_file> --biosystem <biosystem> --cutoff <cutoff_value> --num_cores <num_cores> --out <output_directory>

    Arguments
        --smiles : str, optional
        Input smiles string for prediction. If provided, predictions will be made for this single molecule.

        --query_name : str, optional
        Query name for the provided smiles string. This name will be used in the output file name.

        --file : str, optional
        Input file containing multiple smiles strings and names. The file should be in tab-separated values (TSV) format with columns "name" and "smiles".

        --biosystem : str, optional
        Biosystem to use for predictions. Default is 'all'.

        --cutoff : float, optional
        Cutoff value for predictions. Default is 0.6.

        --num_cores : int, optional
        Number of cores to use for multiprocessing. Default is 4.

        --out : str, optional
        Output directory for predictions. Default is the current directory ('./').

    Description
        The script can be run in two modes:

            - Single Prediction Mode: When the --smiles argument is provided, the script will perform predictions for the single molecule specified by the smiles string.
            - Batch Prediction Mode: When the --file argument is provided, the script will read the input file and perform predictions for each molecule listed in the file using multiprocessing.
    
    Examples
        Single Prediction Mode:
            
            python microberx.py --SMILES "CCO" --query_name "ethanol" --biosystem "all" --cutoff 0.6 --out "./predictions/"
        
            This command will perform predictions for ethanol and save the results in the specified output directory.

        Batch Prediction Mode:
            
            python microberx.py --file "input_molecules.tsv" --biosystem "all" --cutoff 0.6 --num_cores 4 --out "./predictions/"
            
            This command will read the input_molecules.tsv file, perform predictions for each molecule using 4 cores, and save the results in the specified output directory.
            
    Output
        The predictions are saved as TSV files in the specified output directory. Each file is named after the query name provided or the names listed in the input file.      
"""

def main():
    parser = argparse.ArgumentParser(description='MicrobeRX Predictions')
    parser.add_argument('--smiles', type=str, default="", help='Input smiles string for prediction')
    parser.add_argument('--query_name', type=str, default="", help='Query name for the smiles string')
    parser.add_argument('--file', type=str, default="", help='Input file for predictions')
    parser.add_argument('--biosystem', type=str, default='all', help='Biosystem to use for predictions options: all, human, gutmicrobes')
    parser.add_argument('--cutoff', type=float, default=0.6, help='Cutoff for predictions')
    parser.add_argument('--num_cores', type=int, default=4, help='Number of cores to use for multiprocessing')
    parser.add_argument('--out', type=str, default='./', help='Output directory for predictions')

    args = parser.parse_args()

    smiles = args.smiles
    query_name = args.query_name
    file = args.file
    biosystem = args.biosystem
    cutoff = args.cutoff
    num_cores = args.num_cores
    out = args.out

    if smiles:
        query = Chem.MolFromSmiles(smiles)
        
        Predictor = MetabolitePredictor(query, query_name=query_name, biosystem=biosystem, cut_off=cutoff)
        Predictor.run_prediction()
    
        metabolites = Predictor.predicted_metabolites
        
        output_file = os.path.join(out, f"{query_name}.tsv")
        metabolites.to_csv(output_file, sep='\t', index=False)
        print(f"Predictions saved to {output_file}")
    
    elif file:
        mols = pd.read_csv(file, sep="\t")
        
        def _runPredictions(index):
            try:
                name = mols.name[index]
                query = Chem.MolFromSmiles(mols.smiles[index])
                
                Predictor = MetabolitePredictor(query, query_name=name, biosystem=biosystem, cut_off=cutoff)
                Predictor.run_prediction()

                output_file = os.path.join(out, f"{name}.tsv")
                Predictor.predicted_metabolites.to_csv(output_file, sep='\t', index=False)
            
            except Exception as e:
                print(f"Error processing {mols.name[index]}: {e}")
        
        with multiprocessing.Pool(num_cores) as pool:
            list(tqdm(pool.imap(_runPredictions, mols.index), total=len(mols.index)))

if __name__ == '__main__':
    main()