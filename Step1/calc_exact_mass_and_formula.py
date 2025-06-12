import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops
import argparse

# Function to calculate exact mass from SMILES
def calculate_exact_mass(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Descriptors.ExactMolWt(mol)
        else:
            return None
    except:
        return None

# Function to calculate molecular formula from SMILES
def calculate_molecular_formula(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.rdMolDescriptors.CalcMolFormula(mol)
        else:
            return None
    except:
        return None

def main(input_file, output_stem, sample_column, smiles_column):
    # Read the tab-separated input file
    df = pd.read_csv(input_file, sep='\t')
    
    # Check if the specified columns exist
    if sample_column not in df.columns or smiles_column not in df.columns:
        print(f"Error: Specified columns '{sample_column}' or '{smiles_column}' not found in the file.")
        return
    
    # Calculate exact mass and molecular formula for each SMILES in the specified column
    df['ExactMass'] = df[smiles_column].apply(calculate_exact_mass)
    df['MolecularFormula'] = df[smiles_column].apply(calculate_molecular_formula)
    
    # Define output file name
    output_file = f"{output_stem}_exactmass_and_formula.txt"
    
    # Save the updated DataFrame to a new tab-separated file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="This script reads a tab-separated file, identifies sample and SMILES columns, "
                    "calculates the exact mass and molecular formula for each SMILES, and saves the updated data to a new file."
    )
    parser.add_argument("input_file", help="Path to the input tab-separated file.")
    parser.add_argument("output_stem", help="Stem for the output file name. "
                                            "The output file will be named as <output_stem>_exactmass_and_formula.txt.")
    parser.add_argument("--sample_column", default="sample", 
                        help="Name of the column containing sample IDs. Default is 'sample'.")
    parser.add_argument("--smiles_column", default="MOLSMILES", 
                        help="Name of the column containing SMILES strings. Default is 'MOLSMILES'.")
    
    args = parser.parse_args()
    main(args.input_file, args.output_stem, args.sample_column, args.smiles_column)
