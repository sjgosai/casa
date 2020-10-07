import numpy as np
import pandas as pd
import argparse
import math
import csv
from collections import OrderedDict

def get_args():    
    parser = argparse.ArgumentParser(description='Call peaks over CRISPRi screen windows.')
    parser.add_argument('guide_data',help='Input flow-fish guide count data.')
    parser.add_argument('output_file',help='Output file to print TSV to.')
    parser.add_argument('--meta_data','-m',type=str,default=None,help='Tab-delimited metadata table. '+\
                                                                      'Must contain header with `Short_name` '+\
                                                                      'and `cutting_specificity_score` fields.')
    parser.add_argument('--window_size','-w',type=int,default=100,help='Window size over which to subset.')
    parser.add_argument('--subset_frac','-f',type=float,default=0.5,help='Maximum fraction of guides to remove by subsetting.')
    parser.add_argument('--negate','-n',action='store_true',help='Negate `cutting_spcificity_score` to flip sort order.')
    parser.add_argument('--random_seed','-r',type=int,default=None,help='Set random seed for '+\
                                                                        'consistant sampling. '+\
                                                                        'Use current date (e.g., 20200901). '+\
                                                                        'Required if meta_data is None')
    args = parser.parse_args()
    args.step_size = args.window_size
    return args

def check_args(args):
    assert args.subset_frac > 0, "subset_frac aust be greater than 0."
    assert args.subset_frac < 0, "subset_frac aust be less than 1."
    assert args.window_size > 0, "Windows must take up space."
    return True

def check_overlap_bed(interval, array):
    height = array.shape[0]
    intervals = np.stack([np.tile(interval,(height,1)), array],axis=0)
    swaghook = (intervals[0,:,1] < intervals[1,:,1]).astype(int)
    chrom    = (intervals[0,:,0] == intervals[0,:,0])
    overlap  = intervals[1-swaghook,np.arange(height),2] > intervals[swaghook,np.arange(height),1]
    return overlap & chrom

def main(args):
    #####################
    ##   Import data   ##
    #####################
    guide_data = pd.read_table(args.guide_data, sep='\t', header=0, index_col=False)
    if args.meta_data is not None:
        meta_data = pd.read_table(args.meta_data, sep='\t', header=0, index_col='Short_name')
    else:
        np.random.seed(args.random_seed)
        meta_data = {
            'Short_name': guide_data['Coordinates'],
            'cutting_specificity_score': np.random.rand( len(guide_data['Coordinates']) )
        }
        meta_data = pd.DataFrame(meta_data['cutting_specificity_score'], 
                                 index=meta_data['Short_name'])
    if args.negate:
        meta_data['cutting_specificity_score'] = -1 * meta_data['cutting_specificity_score']
    else:
        pass
    ########################################
    ## Fill missing metadata, link guides ##
    ########################################
    spec_list  = []
    for name in guide_data['Coordinates']:
        try:
            spec_list.append( meta_data.loc[name,'cutting_specificity_score'] )
        except KeyError:
            spec_list.append( np.nan )
    guide_data['cutting_specificity_score'] = spec_list
    guide_data = guide_data.fillna( 0.0 )
    ## Split targeting and control guides
    targ_data  = guide_data.loc[ guide_data['Coordinates'].str.contains('chr') ]
    ctrl_data  = guide_data.loc[~guide_data['Coordinates'].str.contains('chr') ]
    ## Parse targeting coordinates
    plus_offsets = [152, 147]
    minus_offsets= [146, 153]
    uniq_chrom= np.unique([coord.split(':')[0] for coord in targ_data['Coordinates']])
    chrom2idx = OrderedDict( [ (x,i) for i,x in enumerate(uniq_chrom) ] )
    idx2chrom = OrderedDict( [ (i,x) for i,x in enumerate(uniq_chrom) ] )
    ## Get targeting positions for each guide
    pos_array = np.array([ ( chrom2idx[coord.split(':')[0]],
                             int(coord.split(':')[1].split('-')[1]) - plus_offsets[0],
                             int(coord.split(':')[1].split('-')[1]) + plus_offsets[1] ) if coord.split(':')[2] == '+' 
                           else ( chrom2idx[coord.split(':')[0]],
                                  int(coord.split(':')[1].split('-')[1]) - minus_offsets[0],
                                  int(coord.split(':')[1].split('-')[1]) + minus_offsets[1] )
                           for coord in targ_data['Coordinates'] ])
    ## Get genomic windows covered on each chrom
    genome_lims = OrderedDict(
        [ (idx, 
           (pos_array[pos_array[:,0] == idx, 1].min(),
            pos_array[pos_array[:,0] == idx, 2].max())
          ) for idx, chrom in idx2chrom.items() ] 
    )
    sliding_window = [ (idx, np.vstack( (np.arange(*lims,args.step_size),
                                        np.minimum(np.arange(*lims,args.step_size)+args.window_size,lims[1])) ).T 
                       )
                       for idx, lims in genome_lims.items() ]
    sliding_window = np.concatenate(
        [ np.concatenate( (np.tile( [[idx]], (a_window.shape[0],1) ), a_window), axis=1 ) 
          for idx, a_window in sliding_window ]
    )
    sliding_window = sliding_window[[ np.any(check_overlap_bed(interval,pos_array)) 
                                      for interval in sliding_window ]]
    ## Get chromosome
    chrom = targ_data['Coordinates'].iloc[0].split(':')[0]
    ## Work through windows and do subset
    guide_filter = []
    for a_window in sliding_window:
        guide_capture = check_overlap_bed(a_window,pos_array)
        current_slice = targ_data.loc[guide_capture]
        capture_count = guide_capture.sum()
        removal_count = math.floor(capture_count*args.subset_frac)
        removal_count = max(0, removal_count - current_slice['Coordinates'].isin(guide_filter).sum())
        current_slice = current_slice.iloc[~current_slice['Coordinates'].isin(guide_filter).values ]
        remove_guides = current_slice.iloc[ current_slice['cutting_specificity_score'].argsort() ] \
                                     .iloc[:removal_count,:] \
                                     .loc[:,'Coordinates'].values
        guide_filter = guide_filter + list(remove_guides)
        
    final_data = pd.concat( [ctrl_data, targ_data.iloc[~targ_data['Coordinates'].isin(guide_filter).values ]] )

    final_data.drop('cutting_specificity_score', axis=1) \
      .to_csv(args.output_file, sep="\t", quoting=csv.QUOTE_NONE, index=False)
    return None

if __name__ == "__main__":
    args = get_args()
    main(args)
