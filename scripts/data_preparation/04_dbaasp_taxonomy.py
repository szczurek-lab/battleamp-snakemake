import argparse
import pandas as pd
import requests
from config import NUM_RECORDS

def fetch_uniprot_taxonomy_info(query, result_dict):
    """
    Fetch taxonomy information from UniProt for a given query and update the result dictionary.

    Parameters:
    - query (str): The query string for taxonomy information.
    - result_dict (dict): Dictionary to store taxonomy information.
    """

    # Construct URL for UniProt taxonomy search
    url = f'https://rest.uniprot.org/taxonomy/search?query={query}&format=tsv'
    # Fetch response from UniProt API
    response = requests.get(url).text
    # Parse response and update result dictionary
    response_dict = dict(
        zip(response.splitlines()[0].split('\t'),
            response.splitlines()[1].split('\t'))
    )
    rank = response_dict['Rank']
    name = response_dict['Scientific name']
    parent = response_dict['Parent']
    # Update result dictionary with taxonomy information
    result_dict[rank] = name
    # Recursively fetch parent taxonomy information
    # If not at the superkingdom level  --> end of taxonomy hierarchy
    if rank != 'superkingdom':
        fetch_uniprot_taxonomy_info(parent, result_dict)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('activity_path')
    parser.add_argument('output_path')
    args = parser.parse_args()
    activity_path = args.activity_path
    output_path = args.output_path

    activity_df = pd.read_csv(activity_path)


    # Filter out targets which were recorded less than NUM_RECORDS times
    activity_df = activity_df.groupby("targetSpecies").filter(lambda x: len(x) > NUM_RECORDS)

    # Initialize a dictionary to store species taxonomy information
    species_taxonomy_dict = {}

    # Iterate through unique target species
    print("Fetching taxonomy information form UniProt")
    for species_entry in activity_df.targetSpecies.unique():
        # Initialize a dictionary to store taxonomy information for each species
        taxonomy_dict = {'dbaasp_entry': species_entry}
        try:
            # First try full name (strain)
            query = species_entry.replace(' ', '+')
            fetch_uniprot_taxonomy_info(query, taxonomy_dict)
        except:
            # If fails (no results in UniProt), try just species name
            # Example: 'Klebsiella pneumoniae KK3 9904' -> 'Klebsiella pneumoniae'
            query = ' '.join(species_entry.split(' ')[:2])
            fetch_uniprot_taxonomy_info(query, taxonomy_dict)
        # Update species taxonomy dictionary with taxonomy information for the current species
        species_taxonomy_dict[species_entry] = taxonomy_dict
    taxonomy_df = pd.DataFrame(species_taxonomy_dict).T
    taxonomy_df = taxonomy_df.rename({'dbaasp_entry': 'strain'}, axis=1)
    taxonomy_df.to_csv(f'{output_path}/dbaasp_taxonomy.csv')