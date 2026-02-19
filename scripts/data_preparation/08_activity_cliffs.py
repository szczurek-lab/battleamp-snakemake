import argparse
import pandas as pd
from pycdhit import cd_hit, read_clstr

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('cluster_file')
    parser.add_argument('activity_file')


    parser.add_argument('output_path')
    args = parser.parse_args()
    cluster_file = args.cluster_file
    activity_file = args.activity_file
    output_path = args.output_path

    activity_df = pd.read_csv(activity_file)
    cluster_df = pd.read_csv(cluster_file)
    cluster_df = cluster_df.drop(['class', 'sequence'], axis=1)
    df_merged = pd.merge(activity_df, cluster_df, left_on='id', right_on='id', how='right')
    print(df_merged.head())

    # In each cluster find a pair of most variable peptides
    unique_sequences = []
    df_list = []

    for cluster in df_merged.cluster.unique():
        cluster_peptides = df_merged[df_merged['cluster'] == cluster].sort_values(by='activity')

        inactive = cluster_peptides[cluster_peptides['class'] == 0]
        active = cluster_peptides[cluster_peptides['class'] == 1]

        if len(inactive) < 1 or len(active) < 1:
            continue

        worst = pd.DataFrame(inactive.sort_values(by='activity').iloc[-1, :]).T
        best = pd.DataFrame(active.sort_values(by='activity').iloc[0, :]).T

        if worst.sequence.item() == best.sequence.item():
            continue

        if worst.sequence.item() in unique_sequences or best.sequence.item() in unique_sequences:
            continue

        unique_sequences.extend([worst.sequence.item(), best.sequence.item()])
        cluster_peptides = pd.concat([best, worst])
        df_list.append(cluster_peptides)

    activity_cliffs = pd.concat(df_list)
    activity_cliffs.to_csv(output_path)