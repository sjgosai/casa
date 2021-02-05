import argparse
import csv
import sys

from collections import OrderedDict

import numpy  as np
import pandas as pd

sys.path.insert(0, "../casa/")
from genome_utils import *
from plot_utils import *

KNOWN_CONSTANTS = {'strandElement': '.', 'guideSpacerSeq': '.', 'guideSeq': '.', 'guideType': '.', 
                   'SeqCounts': 'NaN', 'AnalysisMethod': 'CASA'}

def get_args():
    parser = argparse.ArgumentParser(description='Convert CASA peaks from BED format to the ENCODE standard.')
    parser.add_argument('casa_peak_files', nargs='+', help='CASA output files, raw or post-processed, used as inputs here.')
    parser.add_argument('--casa_guide_file', help='One guide file used as CASA input, if using raw CASA outputs for `casa_peak_files`.')
    parser.add_argument('--gff_file',help='Gene definitions with items in the `geneID` column matching the `gene` column of the CASA post-processed peak files. Unused argument, feature in progress.')
    parser.add_argument('--min_guide_coverage', type=int, default=1, help='Minimum guide coverage for reporting peaks if using raw CASA output files.')
    parser.add_argument('--chrTSS', required=True, help='Chromosome of target transcription start site. (i.e., chrX).')
    parser.add_argument('--startTSS', required=True, help='First promoter nucleotide position.')
    parser.add_argument('--strandGene', required=True, help='Sense strand of target gene.')
    parser.add_argument('--measuredGeneSymbol', required=True, help='Common gene symbol of target gene.')
    parser.add_argument('--measuredEnsemblID', required=True, help='ENSEMBLE ID of target gene.')
    parser.add_argument('--measuredGeneExpression', default='NaN', help='RPKM value or NaN.')
    parser.add_argument('--PerturbationMethod', default='dCas9-KRAB', help='CRISPR based purturbation method.')
    parser.add_argument('--ReadoutMethod', default='FlowFISH', help='Gene expression signal detection method.')
    parser.add_argument('--GenomeBuild', default='hg38', help='Genome build.')
    parser.add_argument('--CellType', required=True, help='Cell-type in which experiment was done.')
    parser.add_argument('--BiosampleGeneticModification', default='NaN', help='Summary of genetic modifications to biosample.')
    parser.add_argument('--BiosampleTreatments', default='NaN', help='Treatment used during experiment.')
    parser.add_argument('--BiosampleTreatmentsAmount', default='NaN', help='Dosage of treatment used.')
    parser.add_argument('--BiosampleTreatmentsDuration', default='NaN', help='Duration of treatment used.')
    parser.add_argument('--output_file', required=True, help='Output file path.')
    return parser.parse_args()

def collect(args):
    test_peaks = [ pd.read_table(
                     fn,sep='\t',header=None,
                     names=['chr','start','end','hdr','pass','strand']
                   ) 
                   for fn in args.casa_peak_files ]
    
    test_peaks = [ merge_bed(
                       df.loc[ df['pass'] == True, ('chr','start','end') ].copy()
                   )
                   for df in test_peaks ]
    
    test_peaks = [ df[ (df['end'] - df['start']) > 100 ]
                   for df in test_peaks ]

    test_peaks = [ df.reset_index(drop=True) for df in test_peaks ]
    
    for i, peak_file in enumerate(test_peaks):
        peak_file['exp_id'] = "{}R{}".format(args.measuredGeneSymbol[:2], i+1)
        peak_file['assay']  = args.measuredGeneSymbol
        peak_file['replicate'] = i+1
    
    test_peaks = pd.concat(test_peaks,axis=0,ignore_index=True).reset_index(drop=True)
    return test_peaks

def merge_and_filter(peaks, args):
    rep_peaks = get_replicating_peaks( peaks )
    fil_peaks = filter_by_guide_coverage( rep_peaks, 
                                          pd.read_table( args.casa_guide_file, sep='\t', header=0 )['Coordinates'],
                                          min_coverage=args.min_guide_coverage
                                        )
    return fil_peaks

def get_bed_format(peaks, args):
    scored_peaks = []
    scores = get_peak_strengths(*args.casa_peak_files)
    for i, a_peak in peaks.iterrows():
        assemble_ = [a_peak['chr'],a_peak['start'],a_peak['end']]
        hits = check_overlap_bed(assemble_, scores.loc[:,('chr','start','end')].values)
        grab_scores = scores.loc[hits, 'score']
        grab_pass   = scores.loc[hits, 'pass']
        grab_sign   = scores.loc[hits, 'sign']
        conflict    = sum([ int(y) if x else 0 for x,y in zip(grab_pass,grab_sign) ]) != grab_pass.sum()
        summit = np.argmax(np.abs(grab_scores))
        peak_score  = grab_scores[summit]
        assemble_.append(peak_score)
        assemble_.append(args.measuredGeneSymbol)
        assemble_.append('.')
        if not conflict:
            scored_peaks.append(assemble_)
    return pd.DataFrame( scored_peaks, columns=['chr','start','end','score','gene','strand'] )

def convert_casa_peaks(peak_df, peak_metadata=None):
    all_keys = ['chrPerturbationTarget', 'startPerturbationTarget',
       'endPerturbationTarget', 'chrTSS', 'startTSS', 'endTSS', 'name',
       'EffectSize', 'strandElement', 'strandGene', 'measuredGeneSymbol',
       'measuredEnsemblID', 'measuredGeneExpression', 'PerturbationTargetID',
       'PerturbationMethod', 'ReadoutMethod', 'AnalysisMethod', 'GenomeBuild',
       'CellType', 'BiosampleGeneticModification', 'BiosampleTreatments',
       'BiosampleTreatmentsAmount', 'BiosampleTreatmentsDuration']
    variable_keys = ['chrPerturbationTarget', 'startPerturbationTarget',
                     'endPerturbationTarget', 'name', 'EffectSize', 'PerturbationTargetID']
    reformatted = []
    for i,x in peak_df.iterrows():
        my_line = {'chrPerturbationTarget': x['chr'], 
                   'startPerturbationTarget': x['start'], 
                   'endPerturbationTarget': x['end'], 
                   'name': '{}:{}-{}:.'.format(x['chr'], x['start'], x['end']), 
                   'EffectSize': x['score'], 
                   'PerturbationTargetID': '{}:{}-{}:.'.format(x['chr'], x['start'], x['end'])}
        all_data= {**my_line, **peak_metadata}
        all_data= OrderedDict([ (my_key, all_data[my_key]) for my_key in all_keys ])
        reformatted.append( pd.Series(all_data) )
    return pd.DataFrame(reformatted)

def main(args):
    with open(args.casa_peak_files[0], 'r') as f:
        line1 = f.readline()
        try:
            _ = float(line1.split()[3])
            RAW_FORMAT = False
        except ValueError:
            line2 = f.readline()
            try:
                _ = float(line2.split()[3])
                RAW_FORMAT = False
            except ValueError:
                RAW_FORMAT = True
        try:
            _ = int(line1.split()[1])
            HEADER_FLAG = None
        except ValueError:
            HEADER_FLAG = 0
    
    if RAW_FORMAT:
        assert args.casa_guide_file is not None, "Must provide guide info. Should have a 'Coordinates' column with entries with format indicating target sequence. E.g., chrX:123456789-123456809:+"
        data = collect(args)
        data = merge_and_filter(data, args)
        data = get_bed_format(data, args)
    else:
        data = pd.read_table(args.casa_peak_files[0], sep='\t', header=HEADER_FLAG)
        data.columns = ['chr','start','end','score','gene','strand']
    
    args.endTSS = args.startTSS
    data = convert_casa_peaks(data, peak_metadata={**KNOWN_CONSTANTS,**vars(args)})
    data.to_csv(args.output_file, sep='\t', index=False, quoting=csv.QUOTE_NONE)
    return data

if __name__ == '__main__':
    args = get_args()
    data = main(args)