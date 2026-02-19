import pandas as pd
import numpy as np

def extract_peptide_data_from_dicts(peptide_id, dict_list):
    # Convert the string representation of the list of dictionaries back to a list of dictionaries
    if isinstance(dict_list, str):
        dict_list = eval(dict_list)  # It got converted to str during saving to .csv

    df = pd.DataFrame()  # Create an empty DataFrame to store activity data

    # Check if the list of dictionaries is not empty
    if dict_list != []:
        try:
            # Convert each dictionary in the list to a DataFrame and store them in a list
            df_list = [pd.DataFrame.from_dict(dict_) for dict_ in dict_list]
        except:
            # If there's an exception during conversion, print the peptide_id and return an empty DataFrame
            print(peptide_id)
            return df

        # Concatenate all DataFrames in df_list into one DataFrame
        df = pd.concat(df_list)

        # Drop the 'description' row if it exists in the index of df
        if 'description' in df.index:
            df = df.drop(['description'], axis=0)

    # Add a column 'id' to the DataFrame to store the peptide_id
    df['id'] = peptide_id

    return df

def write_fasta(headers, sequences, output_file):
    if len(headers) != len(sequences):
        raise ValueError("Number of headers and sequences must be the same.")
    assert len(headers) == len(np.unique(headers))
    with open(output_file, 'w') as fasta_file:
        for header, sequence in zip(headers, sequences):
            fasta_file.write(f">{header}\n{sequence}\n")