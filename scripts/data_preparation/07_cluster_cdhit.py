import argparse
import os
import pandas as pd
from pycdhit import cd_hit, read_clstr
from dbaasp_utils import write_fasta


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('input_fasta')
    parser.add_argument('input_csv')
    parser.add_argument('output_dir')
    parser.add_argument('dataset_id')
    parser.add_argument('min_seq')

    args = parser.parse_args()
    input_fasta = args.input_fasta
    input_csv = args.input_csv
    output_dir = args.output_dir
    dataset_id = args.dataset_id
    min_seq = int(args.min_seq)

    os.makedirs(output_dir, exist_ok=True)

    res = cd_hit(
        i=input_fasta,
        o=f'{output_dir}/{dataset_id}',
        c=0.8,  # clustering threshold
        n=5,  # word size
        sc=1,  # output clusters by decreasing size
    )

    # Convert results to csv
    df_clstr = read_clstr(f"{output_dir}/{dataset_id}.clstr")
    # Take only clusters which have minimum of MIN_NUM_SEQ sequences
    df_clstr['count'] = df_clstr['cluster'].map(
        df_clstr.groupby("cluster")['cluster'].agg(['count'])['count'].to_dict())
    df_clstr = df_clstr[df_clstr['count'] > min_seq]
    df_clstr['id'] = df_clstr['identifier']

    act_df = pd.read_csv(input_csv)
    act_df.index = act_df['id'].astype(str)
    seq_dict = act_df['sequence'].to_dict()
    class_dict = act_df['class'].to_dict()


    # seq_dict.keys = seq.dict()
    df_clstr['sequence'] = df_clstr['id'].map(seq_dict)
    df_clstr['class'] = df_clstr['id'].map(class_dict)
    print('Saving to .csv..')
    df_clstr.to_csv(f'{output_dir}/{dataset_id}_clustered.csv', index=False)