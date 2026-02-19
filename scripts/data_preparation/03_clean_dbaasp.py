import argparse
import pandas as pd
from dbaasp_utils import extract_peptide_data_from_dicts
from config import (
    MIN_LEN,
    MAX_LEN,
    CTERMINUS,
    NTERMINUS,
    STD_AA,
    MIN_CFU,
    MAX_CFU,
    CFU_GROUP,
    ACT_MEASURE,
    MEDIUM,
    PH,
    SALT
)


# Remove sequences with non-standard aa
def check_seq_for_std_aa(test_str):
    try:
        standard_aa = set('ACDEFGHIKLMNPRSTWYQV')
        return set(test_str) <= standard_aa
    except:
        return False


def check_for_terminus_mod(mod_info, allowed_mod):
    mod_info = str(mod_info)
    if mod_info == 'nan':
        return True
    else:
        try:
            mod_info = eval(mod_info)
            mod = mod_info['name']
            return mod in allowed_mod
        except:
            print(mod_info)

def filter_cfu(cfu):
    try:
        cfu = float(cfu)
        return cfu
    except:
        return None

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('dbaasp_path')
    parser.add_argument('output_path')
    args = parser.parse_args()
    dbaasp_path = args.dbaasp_path
    output_path = args.output_path


    dbaasp = pd.read_csv(f'{dbaasp_path}', index_col=0)

    print(f"Number of entries in DBAASP before filtering: {len(dbaasp)}")
    if STD_AA:
        before_filtering = len(dbaasp)
        dbaasp = dbaasp[dbaasp['sequence'].apply(check_seq_for_std_aa)]
        no_dropped = before_filtering - len(dbaasp)
        print(f"Dropped {no_dropped} peptides after removing non-standard amino acids")
    if MIN_LEN:
        before_filtering = len(dbaasp)
        dbaasp = dbaasp[dbaasp['sequence'].apply(lambda x: len(str(x)) >= MIN_LEN)]
        no_dropped = before_filtering - len(dbaasp)
        print(f"Dropped {no_dropped} peptides after removing sequences shorter than {MIN_LEN}")
    if MAX_LEN:
        before_filtering = len(dbaasp)
        dbaasp = dbaasp[dbaasp['sequence'].apply(lambda x: len(str(x)) <= MAX_LEN)]
        no_dropped = before_filtering - len(dbaasp)
        print(f"Dropped {no_dropped} peptides after removing sequences longer than {MAX_LEN}")
    if CTERMINUS:
        before_filtering = len(dbaasp)
        dbaasp = dbaasp[dbaasp['cTerminus'].apply(lambda x: check_for_terminus_mod(x, allowed_mod=CTERMINUS))]
        no_dropped = before_filtering - len(dbaasp)
        print(f"Dropped {no_dropped} peptides after removing sequences with not allowed C-terminus")
    if NTERMINUS:
        before_filtering = len(dbaasp)
        dbaasp = dbaasp[dbaasp['nTerminus'].apply(lambda x: check_for_terminus_mod(x, allowed_mod=NTERMINUS))]
        no_dropped = before_filtering - len(dbaasp)
        print(f"Dropped {no_dropped} peptides after removing sequences with not allowed N-terminus")

    before_filtering = len(dbaasp)
    dbaasp = dbaasp.drop_duplicates(subset=['sequence'], keep=False)
    no_dropped = before_filtering - len(dbaasp)
    print(f"Dropped {no_dropped} peptides after removing duplicate sequences")

    print(f"Number of entries in DBAASP after filtering: {len(dbaasp)}")
    print(f"Saving to {output_path}")
    dbaasp.to_csv(f'{output_path}/dbaasp_sequences.csv')

    # Get activity info per peptide
    print(f"Creating activity dataframe")
    id_list = dbaasp['id'].tolist()
    activity_list = dbaasp['targetActivities'].tolist()

    activity_df = []
    for peptide_id, dict_list in zip(id_list, activity_list):
        activity_df.append(extract_peptide_data_from_dicts(peptide_id, dict_list))
    activity_df = pd.concat(activity_df)
    print(f"Total number of activity measurements: {len(activity_df)}")


    good_cfu_ids = []

    if MIN_CFU:
        min_cfu = activity_df[activity_df['cfu'] != '']
        min_cfu['cfu'] = min_cfu['cfu'].apply(filter_cfu)
        min_cfu = min_cfu[(min_cfu['cfu'] >= float(MIN_CFU))]
        good_cfu_ids.extend(min_cfu.id.tolist())
    if MAX_CFU:
        max_cfu = activity_df[activity_df['cfu'] != '']
        max_cfu['cfu'] = max_cfu['cfu'].apply(filter_cfu)
        max_cfu = max_cfu[(max_cfu['cfu'] <= float(MAX_CFU))]
        good_cfu_ids.extend(max_cfu.id.tolist())
    if CFU_GROUP:
        cfu_group = activity_df[activity_df['cfuGroup'] == CFU_GROUP]
        good_cfu_ids.extend(cfu_group.id.tolist())


    good_cfu_ids = set(good_cfu_ids)
    before_filtering = len(activity_df)
    activity_df = activity_df[activity_df['id'].isin(good_cfu_ids)]
    no_dropped = before_filtering - len(activity_df)
    print(f"Dropped {no_dropped} activity entries after choosing measurements with specified CFU")

    if ACT_MEASURE:
        before_filtering = len(activity_df)
        activity_df = activity_df[activity_df['activityMeasureValue'] == ACT_MEASURE]
        no_dropped = before_filtering - len(activity_df)
        print(f"Dropped {no_dropped} activity entries after choosing measurements with {ACT_MEASURE} as activity measure")

    if MEDIUM:
        before_filtering = len(activity_df)
        activity_df = activity_df[activity_df['medium'].isin(MEDIUM)]
        no_dropped = before_filtering - len(activity_df)
        print(f"Dropped {no_dropped} activity entries after choosing measurements with not allowed medium")

    if PH:
        before_filtering = len(activity_df)
        activity_df = activity_df[activity_df['ph'] == PH]
        no_dropped = before_filtering - len(activity_df)
        print(f"Dropped {no_dropped} activity entries after choosing measurements with specified pH")
    else:
        before_filtering = len(activity_df)
        activity_df = activity_df[activity_df['ph'] == '']
        no_dropped = before_filtering - len(activity_df)
        print(f"Dropped {no_dropped} activity entries after choosing measurements with specified pH")

    if SALT:
        before_filtering = len(activity_df)
        activity_df = activity_df[activity_df['saltType'] == SALT]
        no_dropped = before_filtering - len(activity_df)
        print(f"Dropped {no_dropped} activity entries after choosing measurements with specified salt type")
    else:
        before_filtering = len(activity_df)
        activity_df = activity_df[activity_df['saltType'] == '']
        no_dropped = before_filtering - len(activity_df)
        print(f"Dropped {no_dropped} activity entries after choosing measurements with specified salt type")

    activity_df = activity_df[activity_df['concentration'] != ''] # Concentration must be specified
    activity_df = activity_df.drop([
        'reference',
        'activityMeasureValue',
        'activityMeasureGroup',
        'cfuGroup',
        'ph',
        'saltType',
        'ionicStrength',

    ], axis=1)

    # Total number of clean measurements
    print(f"Total number of clean measurements: {len(activity_df)}")
    activity_df.to_csv(f'{output_path}/dbaasp_activity.csv')