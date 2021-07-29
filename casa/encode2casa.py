import argparse
import csv

from collections import OrderedDict

import numpy  as np
import pandas as pd

HEADER = [
 'chrom',
 'chromStart',
 'chromEnd',
 'name',
 'SeqCounts',
 'strandPerturbationTarget',
 'PerturbationTargetID',
 'chrTSS',
 'startTSS',
 'endTSS',
 'strandGene',
 'measuredGeneSymbol',
 'measuredEnsemblID',
 'guideSpacerSeq',
 'guideSeq',
 'guideType',
 'Notes'
]

def get_args():
    parser = argparse.ArgumentParser(description='Parse ENCODE data for CASA. File must be parsable by pandas.read_csv and contain the following fields: ["chrom", "chromStart", "chromEnd", "strandPerturbationTarget", "SeqCounts", "name"]')
    parser.add_argument('--positive_bin', required=True, help='Guide-wise counts table, in ENCODE standard format, for data collected in the bin where positive signal is enriched.')
    parser.add_argument('--negative_bin', required=True, help='Guide-wise counts table, in ENCODE standard format, for data collected in the bin where negative signal is enriched.')
    parser.add_argument('--accept_dot', action='store_true')
    parser.add_argument('--output', required=True, help='Output file path.')
    return parser.parse_args()

def filter_files(my_file, accept_dot=False):
    if accept_dot:
        hold = my_file.loc[ ~np.isnan(my_file['SeqCounts']) ]
    else:
        hold = my_file.loc[ (my_file['chrom'] != '.') & (~np.isnan(my_file['SeqCounts'])) ]
    return hold

def org_controls(plus_file, minus_file):
    plus_controls = plus_file[ plus_file['guideType'] == 'negative_control' ].copy()
    plus_targeting= plus_file[ plus_file['guideType'] != 'negative_control' ].copy()
    minus_controls = minus_file[ minus_file['guideType'] == 'negative_control' ].copy()
    minus_targeting= minus_file[ minus_file['guideType'] != 'negative_control' ].copy()
    plus_n_nt = ( ~plus_controls['chrom'].astype(str).str.contains('chr') ).sum()
    minus_n_nt = ( ~minus_controls['chrom'].astype(str).str.contains('chr') ).sum()
    print(f'{plus_n_nt} and {minus_n_nt} non-targeting controls found in plus and minus file, rsp.')
    if (plus_n_nt > 100) and (minus_n_nt > 100):
        print('enough NT controls detected, using those')
    else:
        print('not enought NT controls, using ST and STT controls')
        shared_tags = list(
            set(plus_controls['PerturbationTargetID']) \
              .intersection(set(minus_controls['PerturbationTargetID']))
        )
        print(f'Found {len(shared_tags)} useable tags')
        for i, tag in enumerate(list(shared_tags)):
            plus_controls.loc[ plus_controls['PerturbationTargetID'] == tag, 'chrom' ] = 'NaN'
            plus_controls.loc[ plus_controls['PerturbationTargetID'] == tag, 'PerturbationTargetID' ] = f'NT_extra_{i}'
            minus_controls.loc[ minus_controls['PerturbationTargetID'] == tag, 'chrom' ] = 'NaN'
            minus_controls.loc[ minus_controls['PerturbationTargetID'] == tag, 'PerturbationTargetID' ] = f'NT_extra_{i}'
            plus_file = pd.concat([plus_controls, plus_targeting]).reset_index(drop=True)
            minus_file= pd.concat([minus_controls,minus_targeting]).reset_index(drop=True)
    return plus_file, minus_file
        

def combine_pair(plus_file, minus_file):
    Coords_P = [ 
        x['PerturbationTargetID'] if 'chr' not in str(x["chrom"])
        else
        '{}:{}-{}:+'.format(x["chrom"], 
                            int(x["chromStart"])-20, 
                            int(x["chromStart"])
                           ) if '+' in x['strandPerturbationTarget'] 
        else
        '{}:{}-{}:-'.format(x["chrom"], 
                            int(x["chromEnd"])+1, 
                            int(x["chromEnd"])+21
                           ) if '-' in x['strandPerturbationTarget'] 
        else 
        x['PerturbationTargetID']
        for i, x in plus_file.iterrows() 
               ]
    Counts_P = list(plus_file['SeqCounts'].astype('int64'))
    
    Coords_M = [ 
        x['PerturbationTargetID'] if 'chr' not in str(x["chrom"])
        else
        '{}:{}-{}:+'.format(x["chrom"], 
                            int(x["chromStart"])-20, 
                            int(x["chromStart"])
                           ) if '+' in x['strandPerturbationTarget'] 
        else
        '{}:{}-{}:-'.format(x["chrom"], 
                            int(x["chromEnd"])+1, 
                            int(x["chromEnd"])+21
                           ) if '-' in x['strandPerturbationTarget']
        else 
        x['PerturbationTargetID']
        for i, x in minus_file.iterrows() 
               ]
    Counts_M = list(minus_file['SeqCounts'].astype('int64'))
        
    LS_data = pd.DataFrame( zip(Coords_P, Counts_P), columns=['Coordinates','LS_reads'] ).set_index('Coordinates')
    HS_data = pd.DataFrame( zip(Coords_M, Counts_M), columns=['Coordinates','HS_reads'] ).set_index('Coordinates')
    
    pre_join= LS_data.join(HS_data, how='inner')
    merged  = pre_join[~pre_join.index.duplicated()].reset_index(drop=False)
    
    merged.loc[ ~merged['Coordinates'].astype(str).str.contains('chr'), 'Coordinates'] = 'NT'
    
    return merged

def main(args):
    pos_data = pd.read_table(args.positive_bin, sep='\t', header=None, names=HEADER)
    neg_data = pd.read_table(args.negative_bin, sep='\t', header=None, names=HEADER)
    pos_data = filter_files(pos_data, args.accept_dot)
    neg_data = filter_files(neg_data, args.accept_dot)
    pos_data, neg_data = org_controls(pos_data, neg_data)
    merged_data = combine_pair(pos_data, neg_data)
    merged_data.to_csv(args.output, sep='\t', index=False, quoting=csv.QUOTE_NONE)

if __name__ == "__main__":
    args = get_args()
    main(args)
