import argparse
import pandas as pd
import json
import glob

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('jsons_path')
    parser.add_argument('output_path')
    args = parser.parse_args()
    jsons_path = args.jsons_path
    output_path = args.output_path

    dbaasp = pd.DataFrame()
    total = 0
    empty = 0
    print("Reading json files")
    for j_file in glob.glob(f"{jsons_path}/*.json"):
        total += 1
        filename = j_file[j_file.rfind("/") + 1:]
        try:
            with open(j_file) as train_file:
                dict_train = json.load(train_file)
        except json.JSONDecodeError:
            empty += 1
            continue
        peptide = pd.DataFrame.from_dict(dict_train, orient='index')
        dbaasp = pd.concat([dbaasp, peptide], axis=1)
    dbaasp = dbaasp.T
    print("Total number of entries:", total)
    print("Empty entries:", empty)
    print("Total number of valid sequences:", len(dbaasp))

    # Drop entries without sequence
    dbaasp['sequence'] = dbaasp['sequence'].str.strip()
    dbaasp = dbaasp.dropna(subset=['sequence'])
    dbaasp = dbaasp.reset_index().drop('index', axis=1)
    print(f"Saving to {output_path}")
    dbaasp.to_csv(f'{output_path}/dbaasp.csv')