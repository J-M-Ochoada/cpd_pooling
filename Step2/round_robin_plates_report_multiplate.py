import re
import csv
import argparse
import math
from tqdm import tqdm
from collections import defaultdict
from datetime import datetime

# Function to map text-based delimiters to actual delimiter characters
def map_delimiter(delimiter_text):
    if delimiter_text == 'tab':
        return '\t'
    elif delimiter_text == 'space':
        return ' '
    elif delimiter_text == 'comma':
        return ','
    else:
        raise ValueError(f"Unsupported delimiter: {delimiter_text}")

# Function to generate plate and well location
def get_plate_and_well(index, plate_format, timestamp, verbose=False):
    """
    Generate the plate and well location for a given index based on the plate format.
    
    Args:
        index (int): The index of the compound.
        plate_format (int): The format of the plate (96, 384, 1536).
    
    Returns:
        tuple: (plate number, well location as a string).
    """
    if plate_format == 96:
        well_rows = 'ABCDEFGH'  # 8 rows
        num_cols = 12
    elif plate_format == 384:
        well_rows = 'ABCDEFGHIJKLMNOP'  # 16 rows
        num_cols = 24
    elif plate_format == 1536:
        # Generate rows dynamically: A-Z followed by AA-AF
        well_rows = [chr(r) for r in range(ord('A'), ord('Z') + 1)] + \
                    [f"A{chr(r)}" for r in range(ord('A'), ord('F') + 1)]  # A-Z, AA-AF
        num_cols = 48
        
    else:
        raise ValueError("Unsupported plate format. Supported formats are 96, 384, and 1536.")

    well_columns = [f"{i:02}" for i in range(1, num_cols + 1)]
    wells_per_plate = len(well_rows) * len(well_columns)

    plate = (index // wells_per_plate) + 1
    well_index = index % wells_per_plate
    row_index = well_index // len(well_columns)
    col_index = well_index % len(well_columns)

    well_row = well_rows[row_index]  # Access the correct row from the list
    well_column = well_columns[col_index]
    if verbose:
        debug_filename = f"debug_output_plate{plate_format}_{timestamp}.txt"  # Use an f-string for dynamic filename
        with open(debug_filename, "a") as debug_file:
            debug_file.write(f"Index {index}: Row Index {row_index}, Row {well_row}, Column Index {col_index}, Column {well_column}\n")

    return plate, f"{well_row}{well_column}"
    
def read_compounds_from_csv(file_path, sample_column, ExactMass_column, delimiter):
    print(f"Using delimiter: {repr(delimiter)}")  # Debugging: Print the delimiter being used
    compounds = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter=delimiter)
        header = next(reader)  # Read the header

        # Validate and get the indices for the sample and ExactMass columns
        try:
            sample_index = header.index(sample_column)
            ExactMass_index = header.index(ExactMass_column)
        except ValueError:
            raise ValueError(f"CSV file must contain '{sample_column}' and '{ExactMass_column}' columns.")

        # Full header for output
        full_header = header + ["PoolPlate", "PoolWell"]  

        for row in reader:
            # Store the entire row to preserve the order
            compound_data = list(row)
            sample_id = compound_data[sample_index]
            ExactMass = compound_data[ExactMass_index]
            compounds.append((compound_data, sample_id, ExactMass))  # Store the entire row

    return compounds, full_header

def assign_compounds_to_wells(compounds, plate_format, compounds_per_well=None, total_wells=None, verbose=False, timestamp=None):

    # Sort compounds by ExactMass (or another sorting criterion if necessary)
    compounds.sort(key=lambda x: x[2])  # Sort based on sample_id or modify as needed

    total_compounds = len(compounds)
    if total_wells:
        compounds_per_well = math.ceil(total_compounds / total_wells)
    elif compounds_per_well:
        total_wells = math.ceil(total_compounds / compounds_per_well)

    well_assignments = []
    current_well = 0

    # Round-robin assignment of compounds to wells with progress tracking
    for i in tqdm(range(total_compounds), desc="Assigning compounds to wells"):
        # Get the well for the current assignment
        plate, well = get_plate_and_well(current_well, plate_format, timestamp, verbose)

        # Unpack compound data
        compound_data, sample_id, ExactMass = compounds[i]

        # Append the compound data with plate and well, but exclude exact_mass from compound_data
        well_assignments.append((*compound_data, plate, well))  # Add plate and well at the end

        # Move to the next well for the next compound
        current_well = (current_well + 1) % total_wells

    return well_assignments

def calculate_collisions(well_assignments, threshold, sample_column_index, exact_mass_column_index):
    # Group compounds by well and calculate collisions
    well_groups = defaultdict(list)

    for row in well_assignments:
        plate, well = row[-2], row[-1]  # Extract plate and well
        exact_mass = float(row[exact_mass_column_index])  # Extract exact mass based on its column index
        sample_id = row[sample_column_index]  # Extract sample ID based on its column index
        well_key = (plate, well)
        well_groups[well_key].append((exact_mass, sample_id))  # Store (exact_mass, sample_id)

    well_collisions = {}
    comparisons = []  # List to hold exact mass comparisons
    seen_pairs = set()  # To track pairs of comparisons already made

    for well_key, masses in tqdm(well_groups.items(), desc="Calculating collisions"):
        collisions = 0
        for i in range(len(masses)):
            for j in range(i + 1, len(masses)):
                mass_i, id_i = masses[i]
                mass_j, id_j = masses[j]

                # Avoid double counting pairs
                pair = frozenset([id_i, id_j])
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)

                # Compare the masses
                if abs(mass_i - mass_j) <= threshold:
                    collisions += 1
                    # Record the comparison
                    comparisons.append([well_key[0], well_key[1], id_i, mass_i, id_j, mass_j, "Yes"])
                else:
                    comparisons.append([well_key[0], well_key[1], id_i, mass_i, id_j, mass_j, "No"])

        well_collisions[well_key] = collisions

    print(f"Total comparisons made: {len(comparisons)}")  # Debugging: print number of comparisons made
    return well_collisions, comparisons

def output_comparisons_to_csv(comparisons, output_prefix, plate_format, timestamp):
    with open(f'{output_prefix}_exact_mass_comparisons_plate{plate_format}_{timestamp}.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Plate", "Well", "ID_1", "Mass_1", "ID_2", "Mass_2", "Comparison"])

        # Write each comparison in the comparisons list
        for comparison in comparisons:
            writer.writerow(comparison)
    
    print(f"Exact mass comparisons file saved as {output_prefix}_exact_mass_comparisons_plate{plate_format}_{timestamp}.csv")  # Debugging: confirm file is saved

# Function to output the collisions summary to a file
def output_collisions_summary(well_collisions, output_prefix, plate_format, timestamp):
    with open(f'{output_prefix}_plate{plate_format}_collisions_summary_{timestamp}.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Plate", "Well", "Collisions"])

        for well_key, collisions in well_collisions.items():
            plate, well = well_key
            writer.writerow([plate, well, collisions])

def output_wells_to_csv(well_assignments, full_header, output_prefix, plate_format, timestamp):
    # Define the correct row order for sorting
    well_rows = [chr(r) for r in range(ord('A'), ord('Z') + 1)]  # A-Z
    well_rows += [f"A{chr(r)}" for r in range(ord('A'), ord('F') + 1)]  # AA-AF

    # Sort well assignments by Plate first, then Row order, then Column
    def sort_key(row):
        plate = row[-2]  # Plate number
        well = row[-1]  # Well (e.g., A01, AA01, etc.)
        row_part = well[:-2]  # Extract row part (e.g., A, AA)
        col_part = int(well[-2:])  # Extract column part (e.g., 01)
        return (plate, well_rows.index(row_part), col_part)

    well_assignments.sort(key=sort_key)

    # Write the sorted data to the CSV file
    with open(f'{output_prefix}_plate{plate_format}_compounds_{timestamp}.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(full_header)  # Write the full header

        for row in well_assignments:
            writer.writerow(row)  # Write the entire row as it is
def main():
    # Argument parsing and file setup
    parser = argparse.ArgumentParser(description="Assign compounds to plates and wells, and calculate collisions.")
    parser.add_argument('file_path', type=str, help="Path to file containing compounds.")
    parser.add_argument('-c', '--compounds_per_well', type=int, help="Number of compounds per well.")
    parser.add_argument('-w', '--total_wells', type=int, help="Total number of wells.")
    parser.add_argument('-t', '--threshold', type=float, default=0.1, help="Threshold for mass collision detection.")
    parser.add_argument('-o', '--output_prefix', type=str, default='output', help="Prefix for the output CSV files.")
    parser.add_argument('-d', '--delimiter', type=str, default='tab', choices=['tab', 'space', 'comma'], help="Delimiter for the input file (tab, space, comma).")
    parser.add_argument('--sample_column', type=str, default='sample', help="The name of the sample column in the input file.")
    parser.add_argument('--ExactMass_column', type=str, default='ExactMass', help="The name of the ExactMass column in the input file.")
    parser.add_argument('--plate_format', type=int, default=384, choices=[96, 384, 1536], help="Plate format to use (96, 384, or 1536). Default is 384.")
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output.")

    args = parser.parse_args()
    file_path = args.file_path
    compounds_per_well = args.compounds_per_well
    total_wells = args.total_wells
    threshold = args.threshold
    output_prefix = args.output_prefix
    delimiter = map_delimiter(args.delimiter)  # Map the text-based delimiter to actual delimiter
    sample_column = args.sample_column
    ExactMass_column = args.ExactMass_column

    # Generate timestamp
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Read compounds from the CSV file
    compounds, full_header = read_compounds_from_csv(file_path, sample_column, ExactMass_column, delimiter)
    
    # Get the index of the sample and ExactMass columns
    sample_column_index = full_header.index(sample_column)
    exact_mass_column_index = full_header.index(ExactMass_column)

    # Assign compounds to wells
    well_assignments = assign_compounds_to_wells(compounds, args.plate_format, compounds_per_well, total_wells, args.verbose, timestamp)

    # Output wells to CSV
    output_wells_to_csv(well_assignments, full_header, output_prefix, args.plate_format, timestamp)

    # Calculate collisions and output the summary
    well_collisions, comparisons = calculate_collisions(well_assignments, threshold, sample_column_index, exact_mass_column_index)
    output_collisions_summary(well_collisions, output_prefix, args.plate_format, timestamp)

    # Output the exact mass comparisons to a new CSV file
    output_comparisons_to_csv(comparisons, output_prefix, args.plate_format, timestamp)

    # Optional verbose message
    if args.verbose:
        print(f"Script completed successfully. Output files saved with prefix: {output_prefix}, plate format: {args.plate_format}, timestamp: {timestamp}")

if __name__ == "__main__":
    main()
