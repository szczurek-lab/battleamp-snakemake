import argparse
import csv
import os
import io
import gzip
import requests
import pandas as pd

SLAY_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94529/suppl/GSE94529%5FPeplibV3%5FHiseq%5Fproperties%5F020118.csv.gz'

def csv_to_fasta(csv_file, fasta_file):
    with open(csv_file, newline='', encoding='utf-8') as csvfile, open(fasta_file, 'w', encoding='utf-8') as fastafile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            fastafile.write(f">{row['Id']}\n{row['Sequence']}\n")

def main(output_dir):

    os.makedirs(output_dir, exist_ok=True)

    response = requests.get(SLAY_URL)
    response.raise_for_status()

    # Load directly into pandas without saving
    with gzip.open(io.BytesIO(response.content), mode='rt') as f:
        slay_df = pd.read_csv(f)

    # slay_df = pd.read_csv(f'{output_dir}/{GSE_file}', index_col=0)
    slay_df['Sequence'] = slay_df['trunc_peptide']
    slay_df['Id'] = range(len(slay_df))
    slay_df['Id'] = 'slay_' + slay_df['Id'].astype(str)

    slay_df['class'] = 0

    # Peptides with a significant decrease of at least log2 fold −1 were
    # considered to be depleted from the input library and to have potential
    # antimicrobial activity.
    slay_df.loc[slay_df['lfcMLE'] <= -1, 'class'] = 1

    # the 1.7% of the peptide library that did show depletion and
    # potential antimicrobial activity represents 7,968 peptides
    assert sum(slay_df['class'] == 1) == 7968

    slay_df = slay_df[['Id', 'Sequence', 'lfcMLE', 'class']]

    library_len = len(slay_df)
    print(f"Found {library_len} sequences")

    slay_df = slay_df.reset_index(drop=True)
    slay_df = slay_df.dropna(subset='Sequence')
    slay_df = slay_df.drop_duplicates(subset='Sequence')

    print(f'Dropped {library_len - slay_df.shape[0]} sequences (duplicates or empty)')
    slay_df_pos = slay_df[slay_df['lfcMLE'] <= -1].copy()
    slay_df_neg = slay_df[slay_df['lfcMLE'] > -1].copy()

    slay_df.to_csv(f'{output_dir}/slay_all.csv', index=False)
    slay_df_pos.to_csv(f'{output_dir}/slay_positives.csv', index=False)
    slay_df_neg.to_csv(f'{output_dir}/slay_negatives.csv', index=False)

    csv_to_fasta(f'{output_dir}/slay_positives.csv', f'{output_dir}/slay_positives.fasta')
    csv_to_fasta(f'{output_dir}/slay_negatives.csv', f'{output_dir}/slay_negatives.fasta')
    csv_to_fasta(f'{output_dir}/slay_all.csv', f'{output_dir}/slay_all.fasta')

    # Generate labels.tsv for the pipeline
    labels_df = slay_df[['Sequence', 'class']].copy()
    labels_df.columns = ['sequence', 'label']
    labels_df['label'] = labels_df['label'].map({1: 'AMP', 0: 'non-AMP'})
    labels_df.to_csv(f'{output_dir}/labels.tsv', sep='\t', index=False)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('output_dir')
    args = parser.parse_args()
    main(args.output_dir)
