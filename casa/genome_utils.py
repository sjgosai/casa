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
    intervals= np.stack([np.tile(interval,(height,1)), array],axis=0)
    intervals[:,:,1:3] = intervals[:,:,1:3].astype(int)
    swaghook = (intervals[0,:,1] < intervals[1,:,1]).astype(int)
    chrom    = (intervals[0,:,0] == intervals[1,:,0])
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

def get_replicating_peaks(bed_df, use_singletons=False):
    uniq_assays = list(bed_df['assay'].unique())
    assay_reps  = [ bed_df.loc[bed_df['assay'] == assay,'replicate'] \
                      .unique() 
                    for assay in uniq_assays ]
    assay_count = [ len(rep_list) for rep_list in assay_reps ]
    result_peaks= []
    for assay, reps, count in zip(uniq_assays, assay_reps, assay_count):
        if count == 1:
            if use_singletons:
                result_peaks.append( bed_df.loc[ bed_df['assay'] == assay, ('chr','start','end') ] )
            else:
                pass
        else:
            in_assay = []
            assay_sub= bed_df[ bed_df['assay'] == assay ]
            for rep in reps:
                in_rep   = assay_sub[ assay_sub['replicate'] == rep ]
                rep_merge= merge_bed(in_rep)
                in_assay.append( rep_merge )
            assay_merge = pd.concat(in_assay, axis=0).reset_index(drop=True)
            assay_merge = merge_bed( assay_merge, count=True )
            result_peaks.append( assay_merge.loc[ assay_merge['count'] > 1, ('chr','start','end') ] )
    return merge_bed(pd.concat( result_peaks, axis=0 ).reset_index(drop=True))

def extract_txn_starts(gff_df):
    txn_starts_dict = {}
    for i, line in gff_df.iterrows():
        strand = line['strand']
        geneID = line['geneID']
        if strand == '+':
            txn_starts_dict[ geneID ] = (line['chrom'], line['txStart'])
        elif strand == '-':
            txn_starts_dict[ geneID ] = (line['chrom'], line['txEnd'])
    return txn_starts_dict

def write_interact_format(pack_peaks, pack_scores, pack_TSStuple, target_path, ex_tag="HCR_flowFISH"):
    with open(target_path,'w') as f:
        header = ['chrom','chromStart','chromEnd','name','score','value','exp',
                  'color','sourceChrom','sourceStart','sourceEnd', 'sourceName',
                  'sourceStrand','targetChrom','targetStart','targetEnd',
                  'targetName','targetStrand']
        print("\t".join(header),file=f)
        for my_peaks, my_scores, my_TSStuple in zip(pack_peaks, pack_scores, pack_TSStuple):
            tss_id, tss_chr, tss_nt = my_TSStuple
            for a_row in my_peaks.iterrows():
                a_peak = a_row[1]
                full_tag  = "e{}:{}/{}/{}".format(tss_id,a_row[0],tss_id,ex_tag)
                assemble_ = [a_peak['chr'],a_peak['start'],a_peak['end']]
                hits = check_overlap_bed(a_peak.values,
                                         my_scores.loc[:,('chr','start','end')].values)
                grab_scores = my_scores.loc[hits,'score'].values
                summit = np.argmax(np.abs(grab_scores))
                peak_score  = grab_scores[summit]
                assemble_.append(full_tag)               # string name
                assemble_.append(0)                      # uint score
                assemble_.append(peak_score)             # double value
                assemble_.append(ex_tag)                 # string exp
                assemble_.append("#000000")              # string color
                assemble_.append(a_peak['chr'])          # string sourceChrom
                assemble_.append(a_peak['start'])        # uint sourceStart
                assemble_.append(a_peak['end'])          # uint sourceEnd
                assemble_.append(full_tag.split('/')[0]) # string sourceName
                assemble_.append('.')                    # string sourceStrand
                assemble_.append(tss_chr)                # string targetChrom
                assemble_.append(tss_nt)                 # uint targetStart
                assemble_.append(int(tss_nt)+1)          # uint targetEnd
                assemble_.append(full_tag.split('/')[1]) # string targetName
                assemble_.append('.')                    # string targetStrand
                print("\t".join(18*['{}']).format(*assemble_),file=f)
    return None

def write_bed_format(pack_peaks, pack_scores, pack_TSStuple, target_path, ex_tag="HCR_flowFISH"):
    with open(target_path,'w') as f:
        for my_peaks, my_scores, my_TSStuple in zip(pack_peaks, pack_scores, pack_TSStuple):
            tss_id, tss_chr, tss_nt = my_TSStuple
            for a_row in my_peaks.iterrows():
                a_peak = a_row[1]
                full_tag  = "e{}:{}/{}/{}".format(tss_id,a_row[0],tss_id,ex_tag)
                assemble_ = [a_peak['chr'],a_peak['start'],a_peak['end']]
                hits = check_overlap_bed(a_peak.values,
                                         my_scores.loc[:,('chr','start','end')].values)
                grab_scores = my_scores.loc[hits,'score'].values
                summit = np.argmax(np.abs(grab_scores))
                peak_score  = grab_scores[summit]
                assemble_.append(peak_score)
                assemble_.append(full_tag.split('/')[1])
                assemble_.append('.')
                print("\t".join(6*['{}']).format(*assemble_),file=f)
    return None

