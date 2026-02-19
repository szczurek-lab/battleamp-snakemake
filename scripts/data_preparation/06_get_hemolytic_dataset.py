import argparse
import re
import numpy as np
import pandas as pd
from dbaasp_utils import write_fasta, extract_peptide_data_from_dicts
hc50_values = ['50-60% hemolysis', '40-50% hemolysis']


def extract_percentages(series):
    def parse_percentage(x):
        # Prioritize the first number before ±
        std_match = re.search(r'^(\d+\.\d+)\s*(?:±\s*[\d.]+)?%\s*hemolysis', x)
        if std_match:
            return float(std_match.group(1))

        # Then try range format
        range_match = re.search(r'([\d.]+)(?:-([\d.]+))?%\s*hemolysis', x)
        if range_match:
            if range_match.group(2) is None:
                return float(range_match.group(1))
            return (float(range_match.group(1)) + float(range_match.group(2))) / 2

        return None

    return series.apply(parse_percentage)

def parse_peptide(peptide_entry):

    n_measurements = peptide_entry.shape[0]

    # Handle cases with one measurement only
    if n_measurements == 1:
        peptide_entry = peptide_entry.iloc[0, :]
        if peptide_entry['activityMeasureForLysisValue'] == '0% hemolysis':
            return 0
        if peptide_entry['activity'] <= 32 and peptide_entry['percentHemolysis'] > 1:
            return 1
        if peptide_entry['activity'] > 32 and peptide_entry['percentHemolysis'] <= 10:
            return 0
        if peptide_entry['percentHemolysis'] > 50:
            return 1
        if peptide_entry['activityMeasureForLysisGroup'] in hc50_values:  # HC50
            return 1 if peptide_entry['activity'] <= 128 else 0

    # Handle cases with multiple measurements
    else:
        if np.any(peptide_entry['activity'] < 32):
            return 1
        if np.all(peptide_entry['activity'] > 128):
            return 0
        if np.all(peptide_entry['percentHemolysis']) <= 10:
            return 1 if np.all(peptide_entry['activity'] < 32) else 0
        if peptide_entry['activityMeasureForLysisGroup'].isin(hc50_values).any():
            peptide_entry = peptide_entry[peptide_entry['activityMeasureForLysisGroup'].isin(hc50_values)].iloc[0,
                            :]
            return 1 if peptide_entry['activity'] <= 128 else 0
    return None


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('dbaasp_path')
    parser.add_argument('sequence_path')
    parser.add_argument('outdir')
    args = parser.parse_args()
    dbaasp_path = args.dbaasp_path
    sequence_path = args.sequence_path
    outdir = args.outdir

    # Read dfs
    dbaasp = pd.read_csv(dbaasp_path, index_col=0)
    dbaasp_sequences = pd.read_csv(sequence_path, index_col=0)
    dbaasp_sequences = dbaasp_sequences.drop_duplicates(subset='sequence', keep=False)

    dbaasp = dbaasp[dbaasp['id'].isin(dbaasp_sequences['id'].tolist())]
    toxicity_list = dbaasp['hemoliticCytotoxicActivities'].tolist()
    id_list = dbaasp['id'].tolist()

    toxicity_df = []
    for peptide_id, dict_list in zip(id_list, toxicity_list):
        toxicity_df.append(extract_peptide_data_from_dicts(peptide_id, dict_list))
    toxicity_df = pd.concat(toxicity_df)
    print(f"Total number of toxicity measurements: {len(toxicity_df)} for {toxicity_df['id'].nunique()} peptides")

    toxicity_df = toxicity_df[toxicity_df['targetCell'] == 'Human erythrocytes']
    toxicity_df = toxicity_df[toxicity_df['activityMeasureForLysisGroup'] != 'None']
    toxicity_df = toxicity_df[toxicity_df['activity'] != 0]
    print(f"Total number of hemolytic measurements: {len(toxicity_df)} for {toxicity_df['id'].nunique()} peptides")

    toxicity_df['activityMeasureForLysisValue'] = toxicity_df['activityMeasureForLysisValue'].str.lower()
    toxicity_df['activityMeasureForLysisGroup'] = toxicity_df['activityMeasureForLysisGroup'].str.lower()
    toxicity_df['percentHemolysis'] = extract_percentages(toxicity_df['activityMeasureForLysisValue'])

    binary_tox = {}
    for id_ in toxicity_df['id'].unique():
        subset_df = toxicity_df[toxicity_df['id'] == id_]
        binary_tox[id_] = parse_peptide(subset_df)
    toxicity_df['binary_tox'] = toxicity_df['id'].map(binary_tox)

    toxic_ids = toxicity_df[toxicity_df['binary_tox'] == 1]['id'].tolist()
    nontoxic_ids = toxicity_df[toxicity_df['binary_tox'] == 0]['id'].tolist()

    nontoxic_ids = list(set(nontoxic_ids))
    nontoxic_seq = dbaasp[dbaasp['id'].isin(nontoxic_ids)]['sequence'].tolist()

    toxic_ids = list(set(toxic_ids))
    toxic_seq = dbaasp[dbaasp['id'].isin(toxic_ids)]['sequence'].tolist()

    print(f'{len(toxic_ids)} hemolytic peptides')
    print(f'{len(nontoxic_ids)} non-hemolytic peptides')
    fucked = list(set(nontoxic_seq) & set(toxic_seq))

    toxic_outfile = outdir + '/' + f'hemolytic.fasta'
    write_fasta(toxic_ids, toxic_seq, toxic_outfile)

    nontoxic_outfile = outdir + '/' + f'nonhemolytic.fasta'
    write_fasta(nontoxic_ids, nontoxic_seq, nontoxic_outfile)