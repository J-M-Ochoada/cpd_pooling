import argparse
import pandas as pd
import os

def convert_and_split_echo_input(input_file, output_dir, delimiter, dest_plate_col, dest_well_col, sample_col, formula_col):
    # Map user-friendly delimiter names to actual delimiters
    delimiter_map = {
        'comma': ',',
        'space': ' ',
        'tab': '\t'
    }
    # Use the appropriate delimiter based on the argument
    actual_delimiter = delimiter_map.get(delimiter.lower())
    if actual_delimiter is None:
        raise ValueError("Delimiter must be one of: 'comma', 'space', 'tab'.")

    # Load the file with the specified delimiter
    df = pd.read_csv(input_file, delimiter=actual_delimiter)

    # Rename columns to standardized names for easier manipulation
    df = df.rename(columns={
        dest_plate_col: 'Destination_Plate',
        dest_well_col: 'Destination_Well',
        sample_col: 'Name',
        formula_col: 'Formula'
    })

    # Verify required columns are present
    required_columns = ['Destination_Plate', 'Destination_Well', 'Name', 'Formula']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in the input file.")

    # Sort by 'Destination_Plate' and 'Destination_Well' for output order consistency
    df = df.sort_values(by=['Destination_Plate', 'Destination_Well'])

    # Group the data by each unique (Destination_Plate, Destination_Well) combination
    grouped = df.groupby(['Destination_Plate', 'Destination_Well'])

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Header lines to include in each output file
    header_lines = [
        "# All user editable compound information in the PCDL Manager compounds tab can be imported",
        "# A minimum of formula or mass plus one identifier must be populated for compound import (recommended is formula and at least one other identifier, see examples)",
        "# Multiple synonyms can be entered into the Synonyms column when separated by ';' (e.g. L-Isoleucine; 2-Amino-3-methylvaleric acid)",
        "# All compounds will be assumed to be of a neutral ion type, unless a 1 is entered in the Cation or Anion column"
    ]

    # Initialize an index counter
    index = 1  # Start indexing from 1

    # For each group, save a separate file in the specified format
    for (plate, well), group in grouped:
        # Format the filename based on the plate, well, and zero-padded index
        file_name = f"plate_{plate}_well_{well}_{index:04}.csv"  # Zero-pad index to 4 digits
        output_path = os.path.join(output_dir, file_name)

        # Write the header lines and data to the output file
        with open(output_path, 'w') as f:
            for line in header_lines:
                f.write(line + '\n')
            # Write the actual data, selecting only the 'Name' and 'Formula' columns (no 'Mass')
            group[['Name', 'Formula']].to_csv(f, index=False)

        print(f"Saved {output_path}")

        # Increment the index for the next file
        index += 1

        # Reset the index when reaching the maximum wells per plate
        if index > 1536:
            index = 1  # Reset index for the next plate

        # Write the header lines and data to the output file
        with open(output_path, 'w') as f:
            for line in header_lines:
                f.write(line + '\n')
            # Write the actual data, selecting only the 'Name' and 'Formula' columns (no 'Mass')
            group[['Name', 'Formula']].to_csv(f, index=False)

        print(f"Saved {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert and split compound mixtures file into Echo-compatible format per well.")
    parser.add_argument("input_file", help="Path to the input file containing compound mixtures.")
    parser.add_argument("output_dir", help="Directory to save the output files for each plate-well combination.")
    parser.add_argument("-d", "--delimiter", default="comma", choices=["comma", "space", "tab"],
                        help="Delimiter of the input file: 'comma' (default), 'space', or 'tab'.")
    parser.add_argument("--dest_plate", required=True, help="Column name for destination plate.")
    parser.add_argument("--dest_well", required=True, help="Column name for destination well.")
    parser.add_argument("--sample_col", required=True, help="Column name for sample ID.")
    parser.add_argument("--formula_col", required=True, help="Column name for formula.")

    args = parser.parse_args()

    convert_and_split_echo_input(
        args.input_file, 
        args.output_dir, 
        args.delimiter, 
        args.dest_plate, 
        args.dest_well, 
        args.sample_col,
        args.formula_col
    )
