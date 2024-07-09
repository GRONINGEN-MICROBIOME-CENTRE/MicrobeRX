<img src="https://github.com/GRONINGEN-MICROBIOME-CENTRE/MicrobeRX/blob/main/img/logo.png?raw=true"  width="320" height="300">

[**Description**](#description) | [**Requirements**](#requirements) | [**Installation**](#installation) | [**Tutorials**](#tutorials) | [**Workflow**](#workflow) | [**Citation**](#citation) | [**License**](#license) | [**Information**](#information) | [**Disclaimer**](#disclaimer)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8207745.svg)](https://doi.org/10.5281/zenodo.8207745)

![Visitors](https://api.visitorbadge.io/api/visitors?path=https%3A%2F%2Fgithub.com%2FGRONINGEN-MICROBIOME-CENTRE%2FMicrobeRX&label=Views&labelColor=%23697689&countColor=%23ff8a65&style=flat)

- [GitHub](https://github.com/GRONINGEN-MICROBIOME-CENTRE/MicrobeRX)

- [Read the Docs](https://microberx.readthedocs.io/)

## Description

**MicrobeRX is A tool for enzymatic reaction-based metabolite prediction in the gut microbiome.** <br><br>

> The publication and details about the tool can be found here:

[MicrobeRX: A tool for enzymatic-reaction-based metabolite prediction in the gut microbiome]()


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

## Citation

If you use this software or its results in your research, publication, or project, please cite it as follows:

> Ruiz-Moreno AJ, Fu J. MicrobeRX: A tool for enzymatic reaction-based metabolite prediction in the gut microbiome [DOI: 10.5281/zenodo.10204312](https://zenodo.org/record/10204312)

## License
This tool is under GPL-3.0 license, see the LICENSE file for details.

## Information

- This release includes the functionalities to generate reaction rules
- This release includes the functionalities to predict, analyze and visualize metabolites.
- This release includes the functionalities to find microorganism and enzymes.

## Disclaimer 

This software is still under development and may contain bugs or errors. The developers do not guarantee the accuracy, completeness, or reliability of the software or its results. Use it at your own risk and discretion. The software is provided "as is" without any warranty of any kind, either express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and non-infringement. The developers are not liable for any damages, losses, or costs arising from the use of the software or its results.