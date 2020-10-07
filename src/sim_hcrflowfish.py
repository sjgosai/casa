import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

def get_args():    
    parser = argparse.ArgumentParser(description='Simulate HCR Flow-FISH experiment with known knockdown effects.')
    parser.add_argument('--example_cells', required=True, 
                        help='Flow-cytometry readings for cells expressing the HCR target.')
    parser.add_argument('--control_cells', required=True, 
                        help='Flow-cytometry readings for cells with zero expression of the HCR target.')
    parser.add_argument('--plot_path', required=True, 
                        help='Path for data vis plot.')
    parser.add_argument('--output', required=True,
                        help='Output file path')
    parser.add_argument('--target_channel', type=str, default='APC-A', 
                        help='Cytometry channel corresponding to target gene, must be reflected as a header in input files')
    parser.add_argument('--housekeeping_channel', type=str, default='FSC-A', 
                        help='Cytometry channel corresponding to housekeeping gene, must be reflected as a header in input files')
    parser.add_argument('--n_control_guides', type=int, default=1000,
                        help='Number of simulated non-targeting control guides.')
    parser.add_argument('--n_targeting_guides', type=int, default=30,
                        help='Number of simulated guides targeting hypothetical CRE window.')
    parser.add_argument('--n_cre_windows', type=int, default=100,
                        help='Number of simulated active CRE windows.')
    parser.add_argument('--sorting_depth', type=int, default=200,
                        help='Mean number of simulated cells per guide.')
    parser.add_argument('--knockdown_fraction', type=float, default=0.5, 
                        help='Causal guide knockdown effect as a fraction of total possible effect.')
    args = parser.parse_args()
    return args

def check_args(args):
    if os.path.dirname(args.plot_path) != '':
        assert os.path.isdir(os.path.dirname(args.plot_path)), f"Plotting parent directory not found"
    if os.path.dirname(args.output) != '':
        assert os.path.isdir(os.path.dirname(args.output)), f"Output parent directory not found"
    with open(Path(args.example_cells), 'r') as f:
        header = f.readline()
        assert args.target_channel in header, f"Target channel not detected in {args.example_cells} header."
        assert args.housekeeping_channel in header, f"Housekeeping channel not detected in {args.example_cells} header."
    with open(Path(args.control_cells), 'r') as f:
        header = f.readline()
        assert args.target_channel in header, f"Target channel not detected in {args.control_cells} header."
        assert args.housekeeping_channel in header, f"Housekeeping channel not detected in {args.control_cells} header."
    return True

def mvn_mle(data):
    '''data is a (p x n) numpy.array where p is 
    variable dims and n is number of samples.
    
    Returns mean and covariance MLEs
    '''
    assert len(data.shape) == 2, "Only matrices allowed"
    my_mean   = data.mean(1, keepdims=True)
    my_scaled = data - my_mean
    my_cov    = np.matmul( my_scaled, my_scaled.T ) / my_scaled.shape[1]
    return my_mean, my_cov


def simulate_hff(sam_mean, sam_cov, con_mean, con_cov, 
                 target, housekeeping, 
                 n_control_guides=1000, n_targeting_guides=20, 
                 sort_depth=2000, ko_fractions=[0.5]):
    sam_con_delta = sam_mean - con_mean
    control_sets  = np.random.normal(0,max(0,sam_cov[1,1]-con_cov[1,1])**0.5,size=n_control_guides)
    targeting_sets= []
    for i, ko_f in enumerate(ko_fractions):
        targeting_guides= np.random.normal(-ko_f*sam_con_delta[1].item(),
                                           max(0,sam_cov[1,1]-con_cov[1,1])**0.5,size=n_targeting_guides)
        targeting_sets.append(targeting_guides)
    targeting_sets = np.concatenate(targeting_sets)
    
    cell_list = []
    for i, guide in enumerate(control_sets):
        mean_vec = np.array([sam_mean.flatten()[0],sam_mean.flatten()[1]+guide])
        cells = np.exp( np.random.multivariate_normal(mean_vec, con_cov, sort_depth) ).astype(np.int64)
        cells = pd.DataFrame(cells, columns=[housekeeping,target])
        cells['target-norm'] = cells[target] / cells[housekeeping]
        cells['arm'] = 'sim_ctrl'
        cells['ko_frac'] = 0.0
        cells['tag'] = f'NT_{i}'
        cell_list.append(cells.reset_index(drop=True))
        
    for i, guide in enumerate(targeting_sets):
        subset = i // n_targeting_guides
        mean_vec = np.array([sam_mean.flatten()[0],sam_mean.flatten()[1]+guide])
        cells = np.exp( np.random.multivariate_normal(mean_vec, con_cov, sort_depth) ).astype(np.int64)
        cells = pd.DataFrame(cells, columns=[housekeeping,target])
        cells['target-norm'] = cells[target] / cells[housekeeping]
        cells['arm'] = 'sim_targ'
        cells['ko_frac'] = ko_fractions[subset]
        start = 10131+(100000*subset)
        end   = 10131+(100000*subset)+(i%n_targeting_guides)+20
        cells['tag'] = f'chr1:{start}-{end}:+'
        cell_list.append(cells.reset_index(drop=True))
        
    return pd.concat(cell_list, ignore_index=True).reset_index(drop=True)

def simulate_sort(sim_data):
    sim_data['LS_hit'] = sim_data['target-norm'] < np.quantile(sim_data['target-norm'],0.1)
    sim_data['HS_hit'] = sim_data['target-norm'] > np.quantile(sim_data['target-norm'],0.9)
    tag_list = sim_data['tag'].unique()
    low_sort_dict = { tag: sum(sim_data.loc[ sim_data['tag'] == tag, 'LS_hit' ]) 
                      for tag in tag_list }

    high_sort_dict= { tag: sum(sim_data.loc[ sim_data['tag'] == tag, 'HS_hit' ])  
                      for tag in tag_list }

    sort_data = pd.DataFrame(
        [ [low_sort_dict[tag], high_sort_dict[tag]] for tag in tag_list ],
        columns=['LS_reads','HS_reads'],
        index=tag_list
    )
    return sort_data

def plot_data(data, x, y, path, n_data_points=500):
    subsets = []
    for an_arm in data['arm'].unique():
        sub_data = data[ data['arm'] == an_arm ].reset_index(drop=True)
        sub_size = sub_data.shape[0]
        idx_sample = np.random.choice(np.arange(sub_size), size=(n_data_points,), replace=False)
        subsets.append( sub_data.loc[idx_sample] )
    subsets = pd.concat(subsets, ignore_index=True).reset_index(drop=True)
    g = sns.jointplot(data=subsets, x=x, y=y, hue="arm", kind='kde', log_scale=True)
    g.savefig(path)
    return g

def main(args):
    targ_data = pd.read_csv(args.example_cells, header=0)
    ctrl_data = pd.read_csv(args.control_cells, header=0)
    
    targ_data = targ_data[ ~((targ_data[args.target_channel] == 0) | (targ_data[args.housekeeping_channel] == 0)) ]
    ctrl_data = ctrl_data[ ~((ctrl_data[args.target_channel] == 0) | (ctrl_data[args.housekeeping_channel] == 0)) ]

    targ_data['arm'] = 'Sample'
    ctrl_data['arm'] = 'Control'

    all_data = pd.concat([ctrl_data,targ_data],axis=0,ignore_index=True).reset_index(drop=True)
    all_data['target-norm'] = all_data[args.target_channel] / all_data[args.housekeeping_channel]
    
    print("Read input data.", file=sys.stderr)
    
    targ_log = np.log(targ_data.loc[:,(args.housekeeping_channel, args.target_channel)]).values.T
    targ_mean, targ_cov = mvn_mle(targ_log)
    
    ctrl_log = np.log(ctrl_data.loc[:,(args.housekeeping_channel, args.target_channel)]).values.T
    ctrl_mean, ctrl_cov = mvn_mle(ctrl_log)
    
    print("Estimated distribution.", file=sys.stderr)
    
    n_windows = args.n_cre_windows + 1
    sim_cells = simulate_hff(targ_mean, targ_cov, ctrl_mean, ctrl_cov, 
                             args.target_channel, args.housekeeping_channel, 
                             n_control_guides=args.n_control_guides, 
                             n_targeting_guides=args.n_targeting_guides, 
                             sort_depth=args.sorting_depth, 
                             ko_fractions=n_windows*[args.knockdown_fraction])
    
    print("Simulated cells.", file=sys.stderr)
    
    data_vis  = plot_data( pd.concat([all_data, sim_cells]).reset_index(drop=True), 
                           x = args.housekeeping_channel, y = args.target_channel, 
                           path=args.plot_path
                         )
    
    print("Plotted data.", file=sys.stderr)
    
    sim_sort  = simulate_sort(sim_cells)
    
    print("Simulated sorting.", file=sys.stderr)
    
    sim_sort.to_csv(args.output,sep='\t',index_label='Coordinates')
    
    print("Done.", file=sys.stderr)
    
    return sim_sort

if __name__ == '__main__':
    args = get_args()
    check_args(args)
    main(args)










