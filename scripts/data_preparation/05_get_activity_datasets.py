import argparse
import os
import json
import numpy as np
import pandas as pd
from dbaasp_utils import write_fasta
from config import MIN_MIC, MAX_MIC
LEVELS = [None, 'gramplus', 'graminus', 'species', 'strain', 'mic']


def get_seq_by_id(df, id_list):
    selected_df = df.loc[id_list, :]
    seq = selected_df['sequence'].unique().tolist()
    return seq

def check_duplicates(active_sequences, inactive_sequences):
    # Find common sequences
    common_sequences = set(active_sequences) & set(inactive_sequences)
    print(f'Found {len(common_sequences)} sequences that were assigned as both active and inactive')
    assert len(common_sequences) == 0
    print(common_sequences)


def classify_activity(df, min_mic, max_mic, dataset_id):
    print(f'Preparing {dataset_id} binary data')
    df_grouped = df.groupby(['id'])['activity'].agg(['min'])
    active_ids = df_grouped[df_grouped['min'] <= min_mic].index.tolist()
    inactive_ids = df_grouped[df_grouped['min'] >= max_mic].index.tolist()

    df.index = df['id'].tolist()

    active_seq = get_seq_by_id(df, active_ids)
    inactive_seq = get_seq_by_id(df, inactive_ids)


    check_duplicates(active_seq, inactive_seq)
    print(f'{len(active_ids)} active peptides')
    print(f'{len(inactive_ids)} inactive peptides')

    df['class'] = None
    df.loc[active_ids, 'class'] =  1
    df.loc[inactive_ids, 'class'] = 0

    os.makedirs(outdir, exist_ok=True)

    print(f'Saving to csv...')
    df.to_csv(outdir + '/' + f'{dataset_id}.csv', index=False)


    print(f'Saving to fasta...')
    active_outfile = outdir + '/' + f'{dataset_id}_positive.fasta'
    inactive_outfile = outdir + '/' + f'{dataset_id}_negative.fasta'
    all_outfile = outdir + '/' + f'{dataset_id}.fasta'
    df = df.loc[list(set(active_ids) | set(inactive_ids)), ['sequence']]
    df.drop_duplicates(subset=['sequence'], inplace=True)
    all_ids, all_seq = df.index.tolist(), df['sequence'].tolist()
    write_fasta(all_ids, all_seq, all_outfile)
    write_fasta(active_ids, active_seq, active_outfile)
    write_fasta(inactive_ids, inactive_seq, inactive_outfile)

    return

def collect_mic(microbe_df, dataset_id):
    print(f'Preparing {dataset_id} MIC data')

    agg_mic = microbe_df.groupby(['id'])['activity'].agg(['min'])
    agg_mic = agg_mic.reset_index()
    agg_mic.columns = ['id', 'MIC']
    mic_df = pd.merge(agg_mic, microbe_df, left_on='id', right_on='id')

    mic_df = mic_df.drop_duplicates(subset='id')
    mic_df['id'] = 'DBAASP_' + mic_df['id'].astype(str)
    print(f'Found {len(mic_df)} peptides with valid MIC')

    dataset_id = dataset_id + '_mic'
    os.makedirs(outdir, exist_ok=True)

    print(f'Saving to csv...')
    mic_df.to_csv(outdir + '/' + f'{dataset_id}.csv', index=False)
    print('Saving to fasta...')
    headers, sequences = mic_df['id'].tolist(), mic_df['sequence'].tolist()
    outfile = outdir + '/' + f'{dataset_id}.fasta'
    write_fasta(headers, sequences, outfile)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('activity_path')
    parser.add_argument('sequence_path')
    parser.add_argument('taxonomy_path')
    parser.add_argument('moi_path')
    parser.add_argument('outdir')
    args = parser.parse_args()
    activity_path = args.activity_path
    sequence_path = args.sequence_path
    taxonomy_path = args.taxonomy_path
    moi_path = args.moi_path
    _outdir = args.outdir

    # Read dfs
    activity_df = pd.read_csv(f'{activity_path}', index_col=0)
    sequence_df = pd.read_csv(f'{sequence_path}', index_col=0)
    taxonomy_df = pd.read_csv(f'{taxonomy_path}', index_col=0)

    sequence_df = sequence_df.loc[:, [
                               'id',
                               'name',
                               'sequence',
                               'nTerminus',
                               'cTerminus',
                           ]]


    # Combine activity and taxonomic info
    combined_df = pd.merge(taxonomy_df, activity_df, left_on='strain', right_on='targetSpecies')
    combined_df = combined_df.rename({'dbaasp_entry': 'strain'}, axis=1)
    # Clean up
    selected_columns = ['id', 'targetSpecies', 'concentration', 'unit',
                        'medium', 'cfu', 'activity', 'note', 'strain',
                        'species', 'genus', 'family', 'order', 'class',
                        'phylum']
    combined_df = combined_df[selected_columns]
    # Add sequence info
    combined_df = pd.merge(sequence_df[['id', 'name', 'sequence']], combined_df, left_on='id', right_on='id')

    # Read Microbes of Interest (MOI)
    with open(moi_path) as moi_file:
        moi = json.load(moi_file)

    levels = moi.keys()
    for level in levels:
        moi_dict = moi[level]
        target_column = moi_dict['target_column']
        microbes = moi_dict['microbes']
        if level in ['broad', 'gramplus', 'gramminus']:
            outdir = _outdir
            microbe_df = combined_df[combined_df['species'].isin(microbes)]
            dataset_id = level
            classify_activity(
                df=microbe_df,
                min_mic=MIN_MIC,
                max_mic=MAX_MIC,
                dataset_id=dataset_id
            )
        else:
            for microbe in microbes:
                microbe_df = combined_df[combined_df[target_column] == microbe]
                dataset_id = microbe.lower().replace(" ", "")
                outdir =  _outdir + '/' + target_column

                if level == 'mic':
                    collect_mic(microbe_df=microbe_df, dataset_id=dataset_id)
                else:
                    classify_activity(
                        df=microbe_df,
                        min_mic=MIN_MIC,
                        max_mic=MAX_MIC,
                        dataset_id=dataset_id
                    )

