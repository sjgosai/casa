import numpy as np
import pandas as pd

def merge_intervals(intervals, count=False):
    sorted_intervals = intervals[ intervals[:,0].argsort() ]
    merged_intervals = sorted_intervals[0:1]
    counts = [1]
    for i in range(1,sorted_intervals.shape[0]):
        next_interval = sorted_intervals[i:i+1]
        last_interval = merged_intervals[-1:]
        if next_interval[0,0] <= last_interval[0,1]:
            new_max = max( next_interval[0,1], last_interval[0,1] )
            merged_intervals[-1,1] = new_max
            counts[-1] += 1
        else:
            merged_intervals = np.concatenate([merged_intervals,next_interval], 
                                              axis=0)
            counts.append(1)
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
    return pd.concat(hold_, axis=0).reset_index(drop=True)

def check_overlap(interval, array):
    height = array.shape[0]
    intervals = np.stack([np.tile(interval,(height,1)), array],axis=0)
    anchor =  (intervals[0,:,0] < intervals[1,:,0]).astype(int)
    return intervals[1-anchor,np.arange(height),1] > intervals[anchor,np.arange(height),0]

def check_overlap_bed(interval, array):
    height = array.shape[0]
    intervals = np.stack([np.tile(interval,(height,1)), array],axis=0)
    swaghook = (intervals[0,:,1] < intervals[1,:,1]).astype(int)
    chrom    = (intervals[0,:,0] == intervals[0,:,0])
    overlap  = intervals[1-swaghook,np.arange(height),2] > intervals[swaghook,np.arange(height),1]
    return overlap & chrom

def intersect_bed3(array1, array2):
    output = []
    uniq_chr = [ x for x in np.unique(array1[:,0]) if x in np.unique(array2[:,0]) ]
    for on_chr in uniq_chr:
        if sum(array1[:,0] == on_chr) >= sum(array2[:,0] == on_chr):
            longer = array1[ array1[:,0] == on_chr ]
            shorter= array2[ array2[:,0] == on_chr ]
        else:
            longer = array2[ array2[:,0] == on_chr ]
            shorter= array1[ array1[:,0] == on_chr ]
            
        longer =  longer[  longer[:,1].argsort() ]
        shorter= shorter[ shorter[:,1].argsort() ]
        for i in range(shorter.shape[0]):
            left_check = shorter[i,1] < longer[:,2]
            right_check= shorter[i,2] >=longer[:,1]
            overlap = longer[ left_check & right_check ]
            for j in range(overlap.shape[0]):
                left = max(shorter[i,1],overlap[j,1])
                right= min(shorter[i,2],overlap[j,2])
                output.append([shorter[i,0],left,right])
    return pd.DataFrame(output,columns=['chr','start','end'])
