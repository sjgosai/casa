import argparse
import csv

from collections import OrderedDict

import numpy  as np
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description='Parse ENCODE data for CASA. File must be parsable by pandas.read_csv and contain the following fields: ["chrPerturbationTarget", "startPerturbationTarget", "endPerturbationTarget", "strandElement", "SeqCounts", "name"]')
    parser.add_argument('--positive_bin', required=True, help='Guide-wise counts table, in ENCODE standard format, for data collected in the bin where positive signal is enriched.')
    parser.add_argument('--negative_bin', required=True, help='Guide-wise counts table, in ENCODE standard format, for data collected in the bin where negative signal is enriched.')
    parser.add_argument('--output', required=True, help='Output file path.')
    return parser.parse_args()

def combine_pair(plus_file, minus_file):
    Coords_P = [ 
        '{}:{}-{}:+'.format(x["chrPerturbationTarget"], 
                            int(x["startPerturbationTarget"])-20, 
                            int(x["startPerturbationTarget"])
                           ) if '+' in x['strandElement'] 
        else
        '{}:{}-{}:-'.format(x["chrPerturbationTarget"], 
                            int(x["endPerturbationTarget"])+1, 
                            int(x["endPerturbationTarget"])+21
                           ) if '-' in x['strandElement'] 
        else 
        x['name']
        for i, x in plus_file.iterrows() 
               ]
    Counts_P = list(plus_file['SeqCounts'])
    
    Coords_M = [ 
        '{}:{}-{}:+'.format(x["chrPerturbationTarget"], 
                            int(x["startPerturbationTarget"])-20, 
                            int(x["startPerturbationTarget"])
                           ) if '+' in x['strandElement'] 
        else
        '{}:{}-{}:-'.format(x["chrPerturbationTarget"], 
                            int(x["endPerturbationTarget"])+1, 
                            int(x["endPerturbationTarget"])+21
                           ) if '-' in x['strandElement']
        else 
        x['name']
        for i, x in minus_file.iterrows() 
               ]
    Counts_M = list(minus_file['SeqCounts'])
        
    LS_data = pd.DataFrame( zip(Coords_P, Counts_P), columns=['Coordinates','LS_reads'] ).set_index('Coordinates')
    HS_data = pd.DataFrame( zip(Coords_M, Counts_M), columns=['Coordinates','HS_reads'] ).set_index('Coordinates')
    
    merged  = pd.concat([LS_data,HS_data], axis=1, join='inner').reset_index(drop=False)
    
    merged.loc[ ~merged['Coordinates'].str.contains('chr'), 'Coordinates'] = 'NT'
    
    return merged

def main(args):
    pos_data = pd.read_csv(args.positive_bin)
    neg_data = pd.read_csv(args.negative_bin)
    merged_data = combine_pair(pos_data, neg_data)
    merged_data.to_csv(args.output, sep='\t', index=False, quoting=csv.QUOTE_NONE)

if __name__ == "__main__":
    args = get_args()
    main(args)
