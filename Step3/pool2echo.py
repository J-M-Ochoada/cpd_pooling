import argparse
import pandas as pd
import re

def well_to_row_col(well, log_file=None):
    # Match well formats like "A01", "AA01", etc.
    match = re.match(r"([A-Za-z]+)(\d+)", well)
    if match:
        row_letters = match.group(1).upper()  # Row part (e.g., A, AA)
        col_number = int(match.group(2))  # Column part (e.g., 01)
        
        # Convert row letters to a numeric value
        row_number = 0
        for letter in row_letters:
            row_number = row_number * 26 + (ord(letter) - ord('A') + 1)
        
        # Log the conversion to both console and file
        message = f"Converting Well '{well}': Row '{row_letters}' -> {row_number}, Column '{col_number}'"
        print(message)  # Print to console
        if log_file:
            with open(log_file, "a") as f:
                f.write(message + "\n")
        
        return row_number, col_number
    else:
        raise ValueError(f"Invalid well format: {well}")

def convert_to_echo_input(input_file, output_file, delimiter, src_plate, src_well, dest_plate, dest_well, transfer_volume, log_file=None):
    delimiter_map = {
        "comma": ",",
        "tab": "\t",
        "space": " "
    }
    
    if delimiter in delimiter_map:
        delimiter = delimiter_map[delimiter]
    else:
        raise ValueError("Unsupported delimiter. Please use 'comma', 'tab', or 'space'.")

    # Ensure the output file has the .csv extension
    if not output_file.endswith(".csv"):
        output_file = f"{output_file}.csv"

    df = pd.read_csv(input_file, delimiter=delimiter)

    # Rename columns to standardized names
    df = df.rename(columns={
        src_plate: 'Source Plate Name',
        src_well: 'Source Well',
        dest_plate: 'Destination Plate Name',
        dest_well: 'Destination Well'
    })

    # Verify required columns are present
    required_columns = ['Source Plate Name', 'Source Well', 'Destination Plate Name', 'Destination Well']
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in the input file.")

    # Convert Source Well and Destination Well to Row and Column format
    source_row_col = df['Source Well'].apply(lambda well: well_to_row_col(well, log_file))
    destination_row_col = df['Destination Well'].apply(lambda well: well_to_row_col(well, log_file))

    # Concatenate the row and column values into the DataFrame
    df[['Source Row', 'Source Column']] = pd.DataFrame(source_row_col.tolist(), index=df.index)
    df[['Destination Row', 'Destination Column']] = pd.DataFrame(destination_row_col.tolist(), index=df.index)

    # Add the Transfer Volume column
    df['Transfer Volume'] = transfer_volume

    # Reorder columns to match the required output format
    df_output = df[['Source Plate Name', 'Source Column', 'Source Row', 'Destination Plate Name', 
                    'Destination Column', 'Destination Row', 'Transfer Volume']]

    # Save the modified DataFrame to the output file
    df_output.to_csv(output_file, index=False)

    print(f"Converted file saved to {output_file}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a compound mixture file to an Echo-compatible format.")
    parser.add_argument("input_file", help="Path to the input file containing compound mixtures.")
    parser.add_argument("output_file", help="Stem of the output file (the script will add '.csv').")
    parser.add_argument("-d", "--delimiter", choices=["comma", "tab", "space"], default="comma", 
                        help="Delimiter of the input file (choose 'comma', 'tab', or 'space'). Default is 'comma'.")
    parser.add_argument("--src_plate", required=True, help="Column name for source plate.")
    parser.add_argument("--src_well", required=True, help="Column name for source well.")
    parser.add_argument("--dest_plate", required=True, help="Column name for destination plate.")
    parser.add_argument("--dest_well", required=True, help="Column name for destination well.")
    parser.add_argument("--tv", required=True, type=float, help="Transfer volume to be applied to all rows.")
    parser.add_argument("--log_file", help="Path to the log file for well conversions.")

    args = parser.parse_args()

    convert_to_echo_input(
        args.input_file, 
        args.output_file, 
        args.delimiter, 
        args.src_plate, 
        args.src_well, 
        args.dest_plate, 
        args.dest_well,
        args.tv,
        log_file=args.log_file  # Pass the log file
    )