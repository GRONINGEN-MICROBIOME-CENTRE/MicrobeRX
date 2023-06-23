import pandas as pd
from tqdm.notebook import tqdm

def GetStoichiometricSmiles(query,stoichiometry):
    smiles=[]
    name=[]
    if query in compounds.keys():
        i=0
        while i < stoichiometry:
            smiles.append(compounds[query])
            name.append(query)
            i=i+1
    else:pass
    return name,smiles

def GetScheme(string:str=None) -> [str,str]:
    """
    The function takes a chemical equation string as input and returns two strings 
    representing the chemical equation with the coefficients before the chemical species 
    and the second string representing the smiles notation of the compound.
    
    Parameters:
    string: MEtaNetX reaction format.
    """

    equation = re.sub("@\w+", "", string)
    # Split the equation into the left and right sides
    sides = equation.split(" = ")
    left_side, right_side = sides[0], sides[1]

    # Split each side into its individual terms
    left_terms = left_side.split(" + ")
    right_terms = right_side.split(" + ")

    # Iterate through the left side and build the output string
    output = ""
    reaction = ""
    for term in left_terms:
        coefficient, species = term.split(" ")
        smiles=compounds[species]
        output += (species + ".")* int(coefficient)
        reaction+= (smiles + ".")* int(coefficient)
    output = output[:-1]
    reaction = reaction[:-1]#remove last '.'
    output += ">>"
    reaction += ">>"

    # Iterate through the right side and build the output string
    for term in right_terms:
        coefficient, species = term.split(" ")
        smiles=compounds[species]
        output += (species + ".")* int(coefficient)
        reaction+= (smiles + ".")* int(coefficient)
    output = output[:-1]
    reaction = reaction[:-1]#remove last '.'
    
    return output,reaction