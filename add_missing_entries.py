# Lance Tan 6/30/16

import pandas as pd
import re

# Returns tuple ('T_TotalDepth', 'T_RefCount', 'T_AltCount', 'T_AltFreq')
# Example: returns (100, 95, 5, .05) given "DP=100;RD=95;AD=5;VF=.05"
def get_readcounts(readcounts_string):
    pattern = "DP=([\d]+);RD=([\d]+);AD=([\d]+);VF=([.\d]+)"
    match = re.match(pattern=pattern, string=readcounts_string)
    return match.groups()

def main(title_file_path, variants_file_path, outfile):
    print 'Reading title file...'
    # Read title file to get mapping of patient IDs to sample IDs
    title_file = pd.read_table(title_file_path, index_col='Sample_ID')
    patients = title_file.groupby(['Patient_ID']).groups # Dict of {patient ID : [list of sample IDs]}

    normal_samples = title_file[ title_file['Class'].map(lambda sample_class: 'Normal' in sample_class) ]
    normal_sample_ids = set(normal_samples.index)

    del title_file

    print 'Reading variants file...'
    all_calls = pd.read_table(variants_file_path)

    rows_to_add = list()

    for patient_id, sample_ids in patients.iteritems():
        
        # All the variant calls for this patient
        patient_variant_calls = all_calls.loc[all_calls['Sample'].isin(sample_ids)]
        all_call_groups = patient_variant_calls.groupby(['Chrom', 'Start', 'Ref', 'Alt'])
        # Dict of {(chrom, start, ref, alt) : [list of row #s in all_calls]}
        patient_mutations = all_call_groups.groups
        
        # Remove normals from list of samples, if any
        tumor_ids = set(sample_ids) - normal_sample_ids
        
        print 'Looking at patient %s: %d tumor samples, %d mutations' % (patient_id, len(tumor_ids), len(patient_mutations))


        for mutation in patient_mutations.iterkeys():

            mutation_calls = all_call_groups.get_group(mutation)

            # array of sample IDs where this mutation was called
            samples_with_call = mutation_calls.loc[:,'Sample'].values

            # pd.Series of this mutation's DP/RD/AD/VF strings for all of this patient's samples
            readcount_strings = mutation_calls.iloc[0].loc[sample_ids]
            
            row_template = mutation_calls.iloc[0]

            for sample_id in tumor_ids:
                if sample_id not in samples_with_call:
                    # create a new row, change/delete some values, then append
                    new_row = row_template.copy()
                    new_row['Sample'] = sample_id
                    new_row[['FailureReason', 'CallMethod', 'T_Ref+', 'T_Ref-', 'T_Alt+', 'T_Alt-']] = [None] * 6
                    # 4-tuple ('T_TotalDepth', 'T_RefCount', 'T_AltCount', 'T_AltFreq')
                    readcounts = get_readcounts(readcount_strings[sample_id]) 
                    new_row[['T_TotalDepth', 'T_RefCount', 'T_AltCount', 'T_AltFreq']] = readcounts
                    
                    rows_to_add.append(new_row)

    print 'Adding %d new rows...' % len(rows_to_add)
    outdf = all_calls.append(rows_to_add, ignore_index=True)
    print 'Writing to file..'
    outdf.to_csv(outfile, sep='\t', index=False)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Appends entries to annotated_exonic_variants.txt so that for each patient, all of their samples reflects all of their called variants.')
    parser.add_argument(
        "-t",
        "--title_file",
        action="store",
        dest="title_file_path",
        required=True,
        metavar='/path/to/title_file.txt',
        help="Path to title file")
    parser.add_argument(
        "-i",
        "--variants",
        action="store",
        dest="variants_file_path",
        required=True,
        metavar='/path/to/all_annotated_variants.txt',
        help="Path to table of variants")
    parser.add_argument(
        "-o",
        "--out",
        action="store",
        dest="outfile",
        required=True,
        metavar='/path/to/output.txt',
        help="Path to output file")
    args = parser.parse_args()

    main(**vars(args))
