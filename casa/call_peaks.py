import numpy as np
import pandas as pd
import pymc3 as pm, theano.tensor as tt
import os
import sys
import math
import pickle
import argparse
import subprocess
from collections import OrderedDict

def get_args():    
    parser = argparse.ArgumentParser(description='Call peaks over CRISPRi screen windows.')
    parser.add_argument('input_data',help='Input flow-fish count data.')
    parser.add_argument('output_data',help='BED format peak data.')
    parser.add_argument('--job_index','-ji',type=int,default=0,help='Job chunk index. 0 =< job_index < job_range.')
    parser.add_argument('--job_range','-jr',type=int,default=1,help='Number of chunks that peak calling will be split into.')
    parser.add_argument('--window_size','-ws',type=int,default=100,help='Window size for peak calling.')
    parser.add_argument('--step_size','-ss',type=int,default=100,help='Step size for peak calling.')
    parser.add_argument('--rope_threshold','-rt',default=0.693,type=float,help='ROPE threshold for peak calls.')
    parser.add_argument('--no_offsets','-no',action='store_true',help='Use exact coordinates for CRISPR activity. Use if Coordinates provided are exactly the region of effect.')
    args = parser.parse_args()
    return args

def check_args(args):
    assert args.job_range >= 0, "Job range implies no work! Must be greater than 0."
    assert args.job_index < args.job_range, "Job index must be within specified range. Should be in [0, job_range)."
    assert args.window_size > 0, "Windows must take up space."
    assert args.step_size > 0, "Step size must cause window to slide. (E.g. step_size > 0)."
    assert args.step_size <= args.window_size, "Can't have step_size > window_size. Will cause gaps."
    return True

def check_overlap(interval, array):
    height = array.shape[0]
    intervals = np.stack([np.tile(interval,(height,1)), array],axis=0)
    swaghook = (intervals[0,:,0] < intervals[1,:,0]).astype(int)
    return intervals[1-swaghook,np.arange(height),1] > intervals[swaghook,np.arange(height),0]

def check_overlap_bed(interval, array):
    height = array.shape[0]
    intervals = np.stack([np.tile(interval,(height,1)), array],axis=0)
    swaghook = (intervals[0,:,1] < intervals[1,:,1]).astype(int)
    chrom    = (intervals[0,:,0] == intervals[0,:,0])
    overlap  = intervals[1-swaghook,np.arange(height),2] > intervals[swaghook,np.arange(height),1]
    return overlap & chrom

def main(args):
    print("Begin",file=sys.stderr)
    check_args(args)
    #######################################
    ##
    ## Import Data, remove missing guides
    ##
    #######################################
    print("Import Data, remove missing guides",file=sys.stderr)
    data = pd.read_table(args.input_data, sep="\t", header=0)
    hs_zero = data['HS_reads'] > 0
    ls_zero = data['LS_reads'] > 0
    rm_zero = hs_zero & ls_zero
    data = data[ rm_zero ]
    #######################################
    ##
    ## Downsample larger lib to comparible
    ##
    #######################################
    print("Downsample",file=sys.stderr)
    ## Rescale to floats
    rescale = min(data['LS_reads'].sum(),data['HS_reads'].sum()) / data.loc[:,('HS_reads','LS_reads')].sum(axis=0)
    data.loc[:,('HS_reads','LS_reads')] *= rescale
    ## Sample downsized library
    runif              = np.random.uniform(size=data.loc[:,('HS_reads','LS_reads')].shape)
    int_part, sample_p = np.divmod( data.loc[:,('HS_reads','LS_reads')], 1 )
    data.loc[:,('HS_reads','LS_reads')] = int_part + (runif < sample_p)
    ## Return as int
    data.loc[:,('HS_reads','LS_reads')] = data.loc[:,('HS_reads','LS_reads')].astype(int) + 1
    #######################################
    ##
    ## Calc. simple data representations
    ##
    #######################################
    data[ 'beta_mean' ] = data[ 'LS_reads' ] / ( data[ 'LS_reads' ] + data[ 'HS_reads' ] )
    data[ 'log(LS/HS)' ] = np.log( data[ 'LS_reads' ] / data[ 'HS_reads' ] )
    #######################################
    ##
    ## Organize positional information
    ##
    #######################################
    print("Parse positional information",file=sys.stderr)
    ## Line guide effects up to genome
    targ_data = data[ (~data['Coordinates'].str.contains("NT")) &\
                      (~data['Coordinates'].str.contains('CTRL')) &\
                      (~data['Coordinates'].str.contains('FILLER-LV2')) &\
                      (~data['Coordinates'].str.contains('FILLER-SgO')) ]
    if args.no_offsets:
        plus_offsets = [0, 0]
        minus_offsets= [0, 0]
    else:
        plus_offsets = [152, 147]
        minus_offsets= [146, 153]
    uniq_chrom= np.unique([coord.split(':')[0] for coord in targ_data['Coordinates']])
    chrom2idx = OrderedDict( [ (x,i) for i,x in enumerate(uniq_chrom) ] )
    idx2chrom = OrderedDict( [ (i,x) for i,x in enumerate(uniq_chrom) ] )
    pos_array = np.array([ ( chrom2idx[coord.split(':')[0]],
                             int(coord.split(':')[1].split('-')[1]) - plus_offsets[0],
                             int(coord.split(':')[1].split('-')[1]) + plus_offsets[1] ) if coord.split(':')[2] == '+' 
                           else ( chrom2idx[coord.split(':')[0]],
                                  int(coord.split(':')[1].split('-')[1]) - minus_offsets[0],
                                  int(coord.split(':')[1].split('-')[1]) + minus_offsets[1] )
                           for coord in targ_data['Coordinates'] ])
    ## Get genomic windows
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
    #######################################
    ##
    ## Process guide data 
    ##
    #######################################
    print("Process guide data",file=sys.stderr)
    ovl_array = np.stack([ check_overlap_bed(guide_interval,sliding_window) 
                           for guide_interval in pos_array ],axis=0).astype(int)
    ovl_array = np.concatenate((np.zeros_like(ovl_array[:,0:1]),ovl_array),axis=1)
    ovl_dex = pd.DataFrame(ovl_array,columns=["wnd_{}".format(i) for i in np.arange(ovl_array.shape[1])])

    NT_count = data.loc[(data['Coordinates'].str.contains("NT") | data['Coordinates'].str.contains("CTRL")),('Coordinates','HS_reads','LS_reads')].shape[0]
    NT_hold = np.zeros((NT_count,ovl_array.shape[1])).astype(int)
    NT_hold[:,0] = 1
    NT_dex = pd.DataFrame(NT_hold,columns=["wnd_{}".format(i) for i in np.arange(ovl_array.shape[1])])
    
    wind_data = pd.concat((
        pd.concat((data.loc[(data['Coordinates'].str.contains("NT") | data['Coordinates'].str.contains("CTRL")),('Coordinates','HS_reads','LS_reads')].reset_index(drop=True),
               NT_dex.reset_index(drop=True)),axis=1).reset_index(drop=True)
        ,
        pd.concat((targ_data.loc[:,('Coordinates','HS_reads','LS_reads')].reset_index(drop=True),
                   ovl_dex.reset_index(drop=True)),axis=1).reset_index(drop=True)
    ), axis=0, ignore_index=True)
    max_idx = max([ int(item.replace('wnd_','')) for item in wind_data.columns if 'wnd' in item ])
    #######################################
    ##
    ## Call peaks on chunk
    ##
    #######################################
    print("Call peaks",file=sys.stderr)
    chunk_size = math.ceil( float(max_idx) / args.job_range )
    start_idx = 1 + (chunk_size * args.job_index)
    end_idx = start_idx + chunk_size
    
    peak_calls = []
    diff_hdr   = []
        
    for i in range(start_idx,min(max_idx,end_idx)):
        print("Starting wnd_{}".format(i))
        group0 =  (wind_data['wnd_0'] == 1).astype(int)
        group1 =  (wind_data['wnd_{}'.format(i)] == 1).astype(int)
        slicer = np.vstack([group0, group1]).T
        use_data = wind_data[ np.sum(slicer,axis=1) == 1 ]
        slicer = slicer[ np.sum(slicer,axis=1) == 1 ]
        slicer = np.argmax(slicer, axis=1)
        e_mean = np.mean(np.log(wind_data['LS_reads'] / wind_data['HS_reads']))
        e_sd   = np.std(np.log(wind_data['LS_reads'] / wind_data['HS_reads']))
        ct_mean= np.mean(wind_data['LS_reads'].values+ wind_data['HS_reads'].values)
        ct_sd  = np.std(wind_data['LS_reads'].values+ wind_data['HS_reads'].values)
        g_var  = (ct_sd**2) - ct_mean
        if g_var <= 0:
            g_sigma = ct_sd
            print("Warning! Count data is underdispersed, results may be inaccurate.")
        else:
            g_sigma = np.sqrt(g_var)
            
        with pm.Model() as model:
            g = pm.Gamma('guide_intensity',mu=ct_mean,sigma=g_sigma,shape=slicer.shape[0])

            e = pm.Normal('enhancer_activity', mu=e_mean, sigma=e_sd,shape=2)
            p = pm.Deterministic('bin_bias', tt.nnet.sigmoid(e))

            l = pm.Deterministic('low_bin_theta', g*p[slicer] )
            h = pm.Deterministic('high_bin_theta', g*(1-p[slicer]) )

            diff = pm.Deterministic('enhancer_boost', e[1]-e[0])

            l_ct = pm.Poisson('low_reads', mu=l, observed=use_data['LS_reads'])
            h_ct = pm.Poisson('high_reads', mu=h, observed=use_data['HS_reads'])

        with model:
            trace = pm.sample(1000, tune=4000, cores=8)

        hdr = pm.stats.hpd(trace['enhancer_boost'],alpha=0.001)
        thresh   = [-args.rope_threshold,args.rope_threshold]
        the_call = check_overlap(np.array(thresh),np.expand_dims(hdr,axis=0))[0]

        peak_calls.append( the_call )
        diff_hdr.append(hdr)
        
    with open(args.output_data, 'w') as f:
        for i,j in enumerate(range(start_idx,min(max_idx,end_idx))):
            peak_position = sliding_window[j-1]
            region_hdr    = diff_hdr[i]
            region_call   = peak_calls[i] == False
            interval_info = [idx2chrom[peak_position[0]],
                             peak_position[1],peak_position[2],
                             "{},{}".format(*region_hdr),region_call,'.']
            print("{}\t{}\t{}\t{}\t{}\t{}".format(*interval_info),file=f)
        
    print("Done.",file=sys.stderr)
    
if __name__ == "__main__":
    args = get_args()
    main(args)
