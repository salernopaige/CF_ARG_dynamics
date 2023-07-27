import argparse
import os
import pandas as pd

def combine_rgi(sample_files, outdir, columns):
    # Create a dictionary to store combined data for each column pair
    combined_data = pd.DataFrame()
    
    # Iterate over each sample 
    for sample_file in sample_files:
        for column in columns:
            starting_column = column['starting_column']
            renamed_column = column['renamed_column']

            # Extract the common columns
            common_columns = ['ARO Term', 'ARO Accession', 'Reference Model Type', 'Reference DB', 'Resistomes & Variants: Observed Pathogen(s)', 'Drug Class', 'Resistance Mechanism']
            
            try:
                # Check if the file exists and is non-empty
                assert os.path.exists(sample_file), f"Error: File not found - {sample_file}"
                assert os.path.getsize(sample_file) > 0, f"Error: Empty file - {sample_file}"

                # Read the sample file
                sample_data = pd.read_csv(sample_file, delimiter='\t')

                # Get the sample name from the input file name
                sample_name = os.path.splitext(os.path.basename(sample_file))[0].split('.')[0]
                
                # Select the common columns and the renamed column from the current sample table
                selected_columns = common_columns + [starting_column]
                selected_data = sample_data[selected_columns].copy()
                selected_data = selected_data.rename(columns={starting_column: f"{sample_name}_{renamed_column}"})

                # Drop duplicate columns before merging
                final_col_name = f"{sample_name}_{renamed_column}"
                existing_columns = [col for col in selected_data.columns if col not in common_columns]
                duplicate_columns = [col for col in existing_columns if not col.startswith(final_col_name)]
                selected_data = selected_data.drop(columns=duplicate_columns)

            # Error message if processing goes wrong somewhere 
            except Exception as e:
                print(f"Error processing sample file: {sample_file}")
                print(e)

            # Merge the selected data into the combined data for the column pair
            if combined_data.empty:
                combined_data = selected_data
            else:
                combined_data = pd.merge(combined_data, selected_data, on=common_columns, how='outer')
    
    # Save the individual column pair data to a separate file
    column_file = os.path.join(outdir, f'{renamed_column}.csv')
    combined_data.to_csv(column_file, index=False, header=True)
            

if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Combine RGI Output Files to create a new final table')
    parser.add_argument('-i', '--sample_files', nargs='+', help='Paths to sample output files', required=True)
    parser.add_argument('-o', '--outdir', help='Directory for output files', required=True)
    parser.add_argument('-c', '--columns', nargs='+', help='List of starting and renamed columns', required=True)
    args = parser.parse_args()

    # Get the columns from the command-line arguments
    config_columns = []
    for column_pair in args.columns:
        starting_column, renamed_column = column_pair.split(':')
        config_columns.append({'starting_column': starting_column, 'renamed_column': renamed_column})

    # Call the combine_rgi function with the provided arguments
    combine_rgi(args.sample_files, args.outdir, config_columns)
