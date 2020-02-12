import numpy as np
import pandas as pd
import os
import glob
import csv
import argparse

import matplotlib.pyplot as plt
import seaborn as sns

def get_args():    
    parser = argparse.ArgumentParser(description='Summarize CRISPRi screen guide-wise data as a signal track.')
    parser.add_argument('guide_data',help='Input flow-fish guide count data.')
    parser.add_argument('output_file',help='Output file to print track to.')
    parser.add_argument('--summ_plot_tag','-p',type=str,help='File name prefix for data histograms.')
    parser.add_argument('--verbose','-v',action='store_true',help='Print stuff to screen to help with debugging and stuff.')
    parser.add_argument('--median_shift','-m',action='store_true',help='Zero the median of guide-wise scores.')
    args = parser.parse_args()
    return args

def merge_intervals(intervals, count=False):
    sorted_intervals = intervals[ intervals[:,0].argsort() ]
    merged_intervals = sorted_intervals[0:1]
    counts = []
    k = 1
    for i in range(1,sorted_intervals.shape[0]):
        next_interval = sorted_intervals[i:i+1]
        last_interval = merged_intervals[-1:]
        if next_interval[0,0] <= last_interval[0,1]:
            new_max = max( next_interval[0,1], last_interval[0,1] )
            merged_intervals[-1,1] = new_max
            k += 1
        else:
            merged_intervals = np.concatenate([merged_intervals,next_interval], 
                                              axis=0)
            counts.append(k)
            k = 1
    if count:
        return merged_intervals, counts
    else:
        return merged_intervals

def merge_bed(bed_df, count=False):
    bed_ = bed_df.sort_values(['chr','start'],axis=0) \
                 .reset_index(drop=True)
    chr_set = bed_['chr'].unique()
    hold_   = []
    for chrom in chr_set:
        sub_bed = bed_.loc[ bed_['chr'] == chrom, ('start','end') ]
        if count:
            merged_pos, counts = merge_intervals(sub_bed.values, True)
            reorg_df= pd.DataFrame({'chr': chrom, 
                                    'start': merged_pos[:,0],
                                    'end': merged_pos[:,1], 
                                    'count': counts})
        else:
            merged_pos = merge_intervals(sub_bed.values)
            reorg_df= pd.DataFrame({'chr': chrom, 
                                    'start': merged_pos[:,0],
                                    'end': merged_pos[:,1]})
        hold_.append(reorg_df)
    return pd.concat(hold_, axis=0)

def main(args):
    data = pd.read_table(args.guide_data, sep='\t', header=0)
    data.loc[data['Coordinates'].str.contains('NT'), 'Coordinates'] = 'NT'
    data.loc[data['Coordinates'].str.contains('CTRL'),'Coordinates']= 'NT'
    if args.verbose:
        print("Dataframe shape: {}".format(data.shape))
        print("Dataframe top few lines:")
        print(data[0:4])
    hs_zero = data['HS_reads'] > 0
    ls_zero = data['LS_reads'] > 0
    rm_zero = hs_zero & ls_zero
    data = data[ rm_zero ]
    if args.verbose:
        print("Pruned zeros, new shape: {}".format(data.shape))
        print("LS lib size: {}".format(data['LS_reads'].sum()))
        print("HS lib size: {}".format(data['HS_reads'].sum()))
        print("Begin downsample")
    rescale = min(data['LS_reads'].sum(),data['HS_reads'].sum()) / data.loc[:,('HS_reads','LS_reads')].sum(axis=0)
    data.loc[:,('HS_reads','LS_reads')] *= rescale
    ## Sample downsized library
    runif              = np.random.uniform(size=data.loc[:,('HS_reads','LS_reads')].shape)
    int_part, sample_p = np.divmod( data.loc[:,('HS_reads','LS_reads')], 1 )
    data.loc[:,('HS_reads','LS_reads')] = int_part + (runif < sample_p)
    ## Return as int
    data.loc[:,('HS_reads','LS_reads')] = data.loc[:,('HS_reads','LS_reads')].astype(int) + 1
    data['Coordinates2']=data['Coordinates']
    data.loc[data['Coordinates'].str.contains('NT'), 'Coordinates'] = 'NT'
    data[ 'log(LS/HS)' ] = np.log( data[ 'LS_reads' ] / data[ 'HS_reads' ] )
    if args.verbose:
        print("Finished downsample")
        print("LS lib size: {}".format(data['LS_reads'].sum()))
        print("HS lib size: {}".format(data['HS_reads'].sum()))
    if args.summ_plot_tag:
        fig, axes = plt.subplots(2,2,figsize=(15,15))
        for i, tags in enumerate(
            zip(['LS_reads','HS_reads','US_reads'],
                ['Low bin counts','High bin counts', 'Unsorted counts'],
                axes.flat
               )
        ):
            df_key, plot_title, ax = tags
            ax.hist( data.loc[ data[df_key] <= 8000, df_key ], bins=200 )
            ax.set_title(plot_title)
            if df_key == 'US_reads':
                ax.set_xlim(0,2000)
            else:
                ax.set_xlim(0,2000)
        fig.tight_layout()
        fig.savefig("{}__sort_bin_count_hists.pdf".format(args.summ_plot_tag))
        # Activity log-odds
        data_slice_refs = [data[ 'log(LS/HS)' ], 
                           data['log(LS/HS)'][ data['Coordinates'] == 'NT' ],
                           data['log(LS/HS)'][ data['Coordinates'] != 'NT' ]]
        data_slice_tags = ["Marginal", "Non-targeting", "Any-targeting"]
        fig, axes = plt.subplots(3,1,figsize=(10,15))
        for basket in zip(data_slice_refs,data_slice_tags,axes.flat):
            ref, tag, ax = basket
            act_mean = ref.sum()/ref.shape[0]
            act_medi = ref.median()
            act_var  = ref.std()**2
            lead_txt = "{} guide activity distribution. ".format(tag)
            stat_txt = "Mean: {0:0.4f}, Median: {1:0.4f}, Var: {2:0.4f}" \
                         .format(act_mean, act_medi, act_var)
            ax.hist(ref, bins=50)
            ax.set_xlim(-5,5)
            ax.set_title(lead_txt+stat_txt)
        fig.tight_layout()
        fig.savefig("{}__guide_activity_hists.pdf".format(args.summ_plot_tag))
    # Convert targeting data to tracks
    plus_offsets = [152, 147]
    minus_offsets= [146, 153]
    chr_list = [ coord.split(':')[0] for coord in data['Coordinates'] if 'chr' in coord ]
    chr_list = sorted(list(set(chr_list)))
    chr_data = []
    for chrom in chr_list:
        targ_data = data[ data['Coordinates'].str.contains(chrom) ]
        pos_array = np.array([ [ int(coord.split(':')[1].split('-')[1]) - plus_offsets[0],
                                 int(coord.split(':')[1].split('-')[1]) + plus_offsets[1] ] if coord.split(':')[2] == '+' 
                               else [ int(coord.split(':')[1].split('-')[1]) - minus_offsets[0],
                                      int(coord.split(':')[1].split('-')[1]) + minus_offsets[1] ]
                               for coord in targ_data['Coordinates'] ])
        genome_lims = (np.min(pos_array), np.max(pos_array))
        if args.verbose:
            print("{}:{}-{}={}".format(chrom, genome_lims[1], genome_lims[0], genome_lims[1] - genome_lims[0]))
        nt_data = np.zeros( shape=[0, 4] )
        last_start = genome_lims[0]
        last_subset= targ_data[ (pos_array[:,0] <= last_start) & (pos_array[:,1] > last_start) ]
        ep = 1e-200
        merged_pos = merge_intervals(pos_array)
        for i, segment in enumerate( merged_pos ):
            if i == (merged_pos.shape[0]-1):
                segment[1] += 1
            for j, nt_pos in enumerate( range(*segment) ):
                nt_subset = targ_data[ (pos_array[:,0] <= nt_pos) & (pos_array[:,1] > nt_pos) ]
                if not nt_subset['Coordinates'].equals(last_subset['Coordinates']):
                    if np.sum( last_subset[ 'HS_reads' ] ) == 0 or np.sum( last_subset[ 'LS_reads' ] ) == 0:
                        if last_subset.shape[0] >= 1:
                            print(nt_pos)
                            print(last_subset)
                    if last_subset.shape[0] == 1:
                        logit_score = np.log( np.sum(last_subset[ 'LS_reads' ]) / np.sum(last_subset[ 'HS_reads' ]) )
                    elif last_subset.shape[0] > 1:
                        logit_score = np.log( np.sum(last_subset[ 'LS_reads' ]) / np.sum(last_subset[ 'HS_reads' ]) )
                    else:
                        logit_score = np.nan
                    new_entry = np.array([[last_start, last_nt+1, last_subset.shape[0], logit_score]])
                    nt_data   = np.concatenate([nt_data,new_entry])
                    last_start= nt_pos
                last_subset = nt_subset
                last_nt     = nt_pos
            if args.verbose:
                if i % 20 == 0:
                    print("On chunk [{}/{}]".format(i,merged_pos.shape[0]))
        nt_data = pd.DataFrame(
            { 'chr': chrom, 
              'start': nt_data[:,0].astype(np.int64), 
              'end': nt_data[:,1].astype(np.int64), 
              'guide_count': nt_data[:,2].astype(np.int64), 
              'score': nt_data[:,3]
            }
        )
        chr_data.append( nt_data )
    all_data = pd.concat(chr_data, axis=0)
    if args.median_shift:
        norm_factor = data[ 'log(LS/HS)' ].median()
        all_data["score"] = all_data["score"] - norm_factor
    all_data.to_csv(args.output_file, sep="\t", quoting=csv.QUOTE_NONE, index=False)
    return None

if __name__ == "__main__":
    args = get_args()
    main(args)
