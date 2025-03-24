<img src="https://github.com/GRONINGEN-MICROBIOME-CENTRE/MicrobeRX/blob/main/img/logo.png?raw=true"  width="320" height="300">

[**Description**](#description) | [**Requirements**](#requirements) | [**Installation**](#installation) | [**Tutorials**](#tutorials) | [**Workflow**](#workflow) | [**Citation**](#citation) | [**License**](#license) | [**Information**](#information) | [**Disclaimer**](#disclaimer)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8207745.svg)](https://doi.org/10.5281/zenodo.8207745)

![Visitors](https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2FGRONINGEN-MICROBIOME-CENTRE%2FMicrobeRX&label=Views&labelColor=%23697689&countColor=%23ff8a65&style=flat)

- [GitHub](https://github.com/GRONINGEN-MICROBIOME-CENTRE/MicrobeRX)

- [Read the Docs](https://microberx.readthedocs.io/)

## Description

**MicrobeRX is A tool for enzymatic reaction-based metabolite prediction in the gut microbiome.** <br><br>

> The publication and details about the tool can be found here:

[Ruiz-Moreno, A.J., Del Castillo-Izquierdo, Á., Tamargo-Rubio, I. et al. MicrobeRX: a tool for enzymatic-reaction-based metabolite prediction in the gut microbiome. Microbiome 13, 78 (2025).](https://doi.org/10.1186/s40168-025-02070-5)

> Click on the link to watch a video of this project:

[![Video](https://img.youtube.com/vi/RHgsHmAIBkA/0.jpg)](https://www.youtube.com/watch?v=RHgsHmAIBkA)

Question about usage or troubleshooting? Please leave a comment in the discussion section of this repo

## Requirements

MicrobeRX is reliant on a variety of academic software. The following are some of the most noticeable:

- rdkit
- datamol
- pyopenms
- pubchempy
- pandas
- plotly
- mols2grid
- rxnMapper
- ReactionDecoder

## Installation 

**1. Installing MicrobeRX from PIP:**

```
pip install MicrobeRX
```

**2. Installing MicrobeRX from source:**

```
git clone https://github.com/GRONINGEN-MICROBIOME-CENTRE/MicrobeRX.git
cd MicrobeRX
python setup.py install
```

**3. Using MicrobeRX in GoogleColab:**

[Predict with MicrobeRX Colab](https://colab.research.google.com/drive/1bELtC9POifs8ExVqHDoVEu_yaBeWcF8Y?usp=sharing)

[Drug Predictions from MicrobeRX Colab](https://colab.research.google.com/drive/170rIHrZDpXxaL7HjAD-v1LfU5yjADtva?usp=sharing)


## Tutorials

- [Generation of Reaction Rules](https://microberx.readthedocs.io/en/latest/tutorials/ReactionRules.html)
- [Prediction of Metabolites](https://microberx.readthedocs.io/en/latest/tutorials/PredictionMetabolites.html)
- [Omics Integrator](https://microberx.readthedocs.io/en/latest/tutorials/OmicsIntegrator.html)
- [MicrobeRX API](https://microberx.readthedocs.io/en/latest/autoapi/index.html)

## Workflow

<img src="https://github.com/GRONINGEN-MICROBIOME-CENTRE/MicrobeRX/blob/main/img/workflow.png?raw=true"  width="600" height="800">

## Command line

- **Description**
    
    The script can be run in two modes:

        Single Prediction Mode: When the --smiles argument is provided, the script will perform predictions for the single molecule specified by the smiles string.

        Batch Prediction Mode: When the --file argument is provided, the script will read the input file and perform predictions for each molecule listed in the file using multiprocessing.

- **Examples**
        
    Single Prediction Mode:
            ```
            microberx --SMILES "CCO" --query_name "ethanol" --biosystem "all" --cutoff 0.6 --out "./predictions/"
            ```
            
    This command will perform predictions for ethanol and save the results in the specified output directory.
            
    
    
    Batch Prediction Mode:
            ```
            microberx --file "input_molecules.tsv" --biosystem "all" --cutoff 0.6 --num_cores 4 --out "./predictions/"
            ```
            
    This command will read the input_molecules.tsv file, perform predictions for each molecule using 4 cores, and save the results in the specified output directory.


- **Output**

        The predictions are saved as TSV files in the specified output directory. Each file is named after the query name provided or the names listed in the input file.   

- **Arguments**

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
        Output directory for predictions. Default is the current directory ('./')   

## Citation

If you use this software or its results in your research, publication, or project, please cite it as follows:

> [Ruiz-Moreno, A.J., Del Castillo-Izquierdo, Á., Tamargo-Rubio, I. et al. MicrobeRX: a tool for enzymatic-reaction-based metabolite prediction in the gut microbiome. Microbiome 13, 78 (2025).](https://doi.org/10.1186/s40168-025-02070-5)

## License
This tool is under GPL-3.0 license, see the LICENSE file for details.

## Information

- This release includes the functionalities to generate reaction rules
- This release includes the functionalities to predict, analyze and visualize metabolites.
- This release includes the functionalities to find microorganism and enzymes.

## Disclaimer 

This software is still under development and may contain bugs or errors. The developers do not guarantee the accuracy, completeness, or reliability of the software or its results. Use it at your own risk and discretion. The software is provided "as is" without any warranty of any kind, either express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and non-infringement. The developers are not liable for any damages, losses, or costs arising from the use of the software or its results.