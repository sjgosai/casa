import sys

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

from genome_utils import *

#############################
##                         ##
## General use gene drawer ##
##                         ##
#############################

def draw_genes_from_gff(ax, gff, promoter_params=[1000,0.2,1000], gene_ylim=[-7.0,-4.5]):
    """
    Draw genes on your plot.

    Inputs
    ----------
    ax: pyplot.Axes class
        An Axes class that contains the plot you want to draw 
        genes on.
    gff: DataFrame
        A pandas DataFrame containing definitions.
    promoter_params: list
        Three parameters that affect how promoters are drawn for 
        each gene. First is the horizontal length of the promoter 
        arrow in number of nucleotides. Second is arrow width in 
        terms of the y-axis varaible being plotted. Third is arrow 
        head length in number of nucleotides.

    Returns
    -------
    None
    
    """
    assert gene_ylim[0] <= gene_ylim[1], "Provide gene_ylim as [lower_val, higher_val]"
    gene_block = [gene_ylim[0], (0.8 * gene_ylim[1]) + (0.2 * gene_ylim[0])]
    prom_height= [gene_block[1], gene_ylim[1]]
    for i, line in gff.iterrows():
        exonStarts = [int(x) for x in line['exonStarts'].split(',')]
        exonEnds   = [int(x) for x in line['exonEnds'].split(',')]
        direction  = 2*int(line['strand'] == '+') - 1
        for start, end in zip(exonStarts, exonEnds):
            ax.fill_betweenx(gene_block,[start,start],[end,end],facecolor='black')
        ax.hlines(y=np.mean(gene_block),xmin=exonStarts[0],xmax=exonEnds[-1])
        if direction == 1:
            x_loc = exonStarts[0]
        else:
            x_loc = exonEnds[-1]
        ax.vlines(x=x_loc,ymin=prom_height[0],ymax=prom_height[1])
        ax.arrow(x=x_loc,y=prom_height[1],dx=direction*promoter_params[0],dy=0, 
                 head_width=promoter_params[1], head_length=promoter_params[2], fc='k', ec='k')
    return None

################################
##                            ##
## Guide and Peak combo plots ##
##                            ##
################################

def plot_hff_cutsites(plot_interval, cutsite_data, peak_data, plot_ids, get_chrom=None):
    """
    Plot guide-wise assay summaries for individual replicates at the 
    guide cut-sites. Include a track of significant peaks above guide 
    score plot. Important: Before using this script, cutsite_data and 
    peak_data should be filtered to a single chromosome.

    Inputs
    ----------
    plot_interval: list
        A two integer list that defines the genomic field of view.
    cutsite_data: DataFrame
        A pandas DataFrame containing guide-wise activity scores for 
        the various replicates to be tested.
    peak_data: DataFrame
        A pandas DataFrame similar to a BED format file. This specifies 
        DataFrame specfies the significant peaks that are called in each 
        replicate.
    plot_ids:
        A list of IDs that correspond to column names in cutsite_data and 
        peak_data. These IDs will be used to extract relevent data.

    Returns
    -------
    (fig, ax): tuple of pyplot objects
        A tuple containing a pyplot.Figure class and a pyplot.Axis class.
    """
    if get_chrom is not None:
        cutsite_data = cutsite_data[ cutsite_data.index.str.contains(get_chrom+":") ]
        peak_data    = peak_data[ peak_data['chr'] == get_chrom ]
    # Subset cutsite scores
    plot_id_slicer = [an_id for an_id in plot_ids if an_id in cutsite_data.columns]
    sub_cuts = cutsite_data.loc[:,plot_id_slicer]
    sub_cuts['cutsite'] = [ int(coord.split(':')[1].split('-')[1]) - 4 
                            if coord.split(':')[-1] == '+' 
                            else 
                            int(coord.split(':')[1].split('-')[0]) + 3 
                            for coord in sub_cuts.index ]
    slice_cuts = check_overlap(plot_interval,np.vstack([cutsite_data['cutsite'].values, 
                                                       (cutsite_data['cutsite']+1).values]).T)
    sub_cuts = sub_cuts.loc[slice_cuts, :]
    # Subset peak intervals
    sub_peaks= check_overlap(plot_interval, peak_data.loc[:,('start','end')].values)
    sub_peaks= peak_data.loc[sub_peaks,:]
    sub_peaks= sub_peaks.loc[ sub_peaks['exp_id'].isin(plot_ids) ]
    
    cut_types  = np.unique(sub_cuts.columns)
    peak_types = np.unique(sub_peaks['exp_id'])
    
    color_cycle= plt.rcParams['axes.prop_cycle'].by_key()['color']
    cycle_len  = len(color_cycle)
    col_dict = { exp_id: color_cycle[idx % cycle_len] 
                 for idx, exp_id
                 in enumerate(plot_ids) }
    
    score_max = np.nanmax(sub_cuts.loc[:, plot_id_slicer].values)
    score_min = np.nanmin(sub_cuts.loc[:, plot_id_slicer].values)
    score_gap = score_max - \
                score_min
    fig = plt.figure(figsize=(12,6))
    ax  = plt.subplot(111)
    for i, exp_id in enumerate(col_dict.keys()):
        sub_sub_peaks = sub_peaks.loc[ sub_peaks['exp_id'] == exp_id, : ]
        peak_position = score_max + ( ( 0.2+(i*0.05) ) * score_gap )
        for j, row in sub_sub_peaks.iterrows():
            ax.hlines(y=peak_position, 
                      xmin=row['start'], xmax=row['end'],
                      color=col_dict[exp_id])
        if exp_id in cut_types:
            null_filter = sub_cuts[exp_id].isnull()
            ax.scatter(sub_cuts.loc[~null_filter,'cutsite'].values,
                       sub_cuts.loc[~null_filter,exp_id].values,
                       color=col_dict[exp_id],s=8,alpha=0.5)
            
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(*plot_interval)
    custom_lines = [ Line2D([0], [0], color=col_dict[color]) for color in plot_ids ]
    ax.legend(custom_lines, plot_ids,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    return fig, ax

def plot_combined_cutsites(plot_interval, cutsite_data, peak_data, plot_ids, merge_style='replicate', get_chrom=None):
    if get_chrom is not None:
        cutsite_data = cutsite_data[ cutsite_data.index.str.contains(get_chrom+":") ]
        peak_data    = peak_data[ peak_data['chr'] == get_chrom ]
    # Subset cutsite scores
    plot_id_slicer = [an_id for an_id in plot_ids if an_id in cutsite_data.columns]
    sub_cuts = cutsite_data.loc[:,plot_id_slicer]
    sub_cuts['cutsite'] = [ int(coord.split(':')[1].split('-')[1]) - 4 
                            if coord.split(':')[-1] == '+' 
                            else 
                            int(coord.split(':')[1].split('-')[0]) + 3 
                            for coord in sub_cuts.index ]
    slice_cuts = check_overlap(plot_interval,np.vstack([cutsite_data['cutsite'].values, 
                                                       (cutsite_data['cutsite']+1).values]).T)
    sub_cuts = sub_cuts.loc[slice_cuts, :]
    # Subset peak intervals
    sub_peaks= check_overlap(plot_interval, peak_data.loc[:,('start','end')].values)
    sub_peaks= peak_data.loc[sub_peaks,:]
    sub_peaks= sub_peaks.loc[ sub_peaks['exp_id'].isin(plot_ids) ]
    
    cut_types  = np.unique(sub_cuts.columns)
    peak_types = np.unique(sub_peaks['exp_id'])
    
    score_max = np.nanmax(sub_cuts.loc[:, plot_id_slicer].values)
    score_min = np.nanmin(sub_cuts.loc[:, plot_id_slicer].values)
    score_gap = score_max - \
                score_min
    fig = plt.figure(figsize=(12,6))
    ax  = plt.subplot(111)
    avail_data = peak_data.loc[peak_data['exp_id'].isin(plot_ids),('exp_id','assay','replicate')].drop_duplicates()
    exp2assay = {}
    assay2exp = {}
    for row in avail_data.iterrows():
        exp2assay[row[1]['exp_id']] = row[1]['assay']
        try:
            assay2exp[row[1]['assay']].append(row[1]['exp_id'])
        except:
            assay2exp[row[1]['assay']] = [row[1]['exp_id']]
    for i, assay in enumerate(avail_data['assay'].unique()):
        color = plt.rcParams['axes.prop_cycle'].by_key()['color'][i]
        peak_position = score_max + ( ( 0.2+(i*0.05) ) * score_gap )
        if merge_style == 'overlap':
            for exp_id in assay2exp[assay]:
                sub_sub_peaks = sub_peaks.loc[ sub_peaks['exp_id'] == exp_id, : ]
                for j, row in sub_sub_peaks.iterrows():
                    ax.hlines(y=peak_position, 
                              xmin=row['start'], xmax=row['end'],
                              color=color)
        elif merge_style == 'replicate':
            collect_merge = []
            for exp_id in assay2exp[assay]:
                exp_data  = sub_peaks.loc[ sub_peaks['exp_id'] == exp_id, : ]
                if assay == 'GATA1':
                    print(exp_data)
                if exp_data.shape[0] > 0:
                    sub_merge = merge_bed(exp_data)
                    collect_merge.append( sub_merge )
            merge_reps = pd.DataFrame()
            if len(collect_merge)  == 1:
                merge_reps = collect_merge[0]
            elif len(collect_merge) > 1:
                collect_merge = pd.concat(collect_merge,axis=0).reset_index(drop=True)
                merge_reps = merge_bed(collect_merge, count=True)
                merge_reps = merge_reps[ merge_reps['count'] > 1 ]
            for j, row in merge_reps.iterrows():
                ax.hlines(y=peak_position, xmin=row['start'], xmax=row['end'], color=color)
        else:
            print("Peak merging style {} not implmented".format(merge_style))
            raise ValueError
        have_cuts = [exp_id for exp_id in assay2exp[assay] if exp_id in sub_cuts.columns]
        if len(have_cuts) > 0:
            print(assay)
            ax.scatter(sub_cuts['cutsite'],sub_cuts.loc[:,have_cuts].mean(axis=1),color=color,s=4,alpha=0.4)
            ax.vlines(x=sub_cuts['cutsite'],
                      ymin=sub_cuts.loc[:,have_cuts].min(axis=1),
                      ymax=sub_cuts.loc[:,have_cuts].max(axis=1),
                      color=color,linewidth=0.5,alpha=0.2)
            print(sum((sub_cuts.loc[:,have_cuts].max(axis=1) - sub_cuts.loc[:,have_cuts].min(axis=1)) < 0.001))
            
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    custom_lines = [ Line2D([0], [0], color=plt.rcParams['axes.prop_cycle'].by_key()['color'][i]) for i, assay in enumerate(avail_data['assay'].unique()) ]
    ax.legend(custom_lines, avail_data['assay'].unique(),bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.set_xlim(*plot_interval)
    plt.tight_layout()
    return fig, ax

#####################
##                 ##
## Connect-o-grams ##
##                 ##
#####################

def gff_to_locus_lims(gff, rel_pad=0.0):
    hard_min = min([ int(j.split(',')[0]) for j in gff['exonStarts'] ])
    hard_max = max([ int(j.split(',')[-1]) for j in gff['exonEnds'] ])
    gap = hard_max - hard_min
    padding = int(rel_pad * gap)
    return [ hard_min - padding, hard_max + padding ]

def draw_bed_blocks(ax, bed, ylims=[1.25,1.75]):
    for i, line in bed.iterrows():
        start = int(line['start'])
        end   = int(line['end'])
        ax.fill_betweenx(ylims, start, end, facecolor='black')
    return None

def connect_bed_to_genes(ax, bed, gene_target, y_anchor=1.25, y_target=1.0, score_bed=None, xlims=None):
    for i, line in bed.iterrows():
        start = int(line['start'])
        end   = int(line['end'])

        if score_bed is None:
            alph = 1.0
            col  = 'black'
        else:
            midpt = score_bed['score'].median()
            cap   = max(np.abs(score_bed['score'] - midpt))
            slicer= check_overlap([start,end],score_bed.loc[:,('start','end')].values)
            ovl   = score_bed[ slicer ]
            best  = ovl.iloc[ (ovl['score'] - midpt).abs().values.argmax() ]['score']
            if ovl['pass'].sum() > 0:
                alph = np.abs(best-midpt) / cap
                if best >= midpt:
                    col = 'black'
                else:
                    col = 'red'
            else:
                alph = 0.0
                col  = 'black'
        
        block_anchor = start + int( (end-start)/2 )
        dx = gene_target - block_anchor
        dy = y_target - y_anchor
        
        if xlims is None:
            DRAW_FLAG = True
        elif xlims[0] <= block_anchor and xlims[1] >= block_anchor and \
             xlims[0] <= gene_target and xlims[1] >= gene_target:
            DRAW_FLAG = True
        else:
            DRAW_FLAG = False
        
        if DRAW_FLAG:
            ax.annotate('',xy=(gene_target,y_target),xytext=(block_anchor,y_anchor),
                        arrowprops=dict(arrowstyle="->",alpha=alph,color=col))
    return None