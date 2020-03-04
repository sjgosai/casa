import os
import sys
import argparse
import subprocess
import tempfile

def get_args():
    parser = argparse.ArgumentParser(description='Call peaks over CRISPRi screen windows.')
    parser.add_argument('input_data',help='Input flow-fish count data.')
    parser.add_argument('output_tag',help='BED format peak data.')
    parser.add_argument('--billed_project','-b',type=str,required=True,help='GCP project ID for billing.')
    parser.add_argument('--gs_bucket','-g',type=str,required=True,help='Valid Google Storage Bucket path for which you have read/write access.')
    parser.add_argument('--window_size','-ws',type=int,default=100,help='Window size for peak calling.')
    parser.add_argument('--step_size','-ss',type=int,default=100,help='Step size for peak calling.')
    parser.add_argument('--rope_threshold','-rt',default=0.693,type=float,help='ROPE threshold for peak calls.')
    parser.add_argument('--no_offsets','-no',action='store_true',help='Use exact coordinates for CRISPR activity. Use if Coordinates provided are exactly the region of effect.')
    parser.add_argument('--job_count','-j',type=int,default=1,help='Number of jobs to split analysis over.')
    parser.add_argument('--compute_zones','-z',type=str,default='us-*',help='Zones for VMs. Remember, costs very between zones.')
    parser.add_argument('--preemptable','-p',action='store_true',help='Flag to use preemptable VMs.')
    args = parser.parse_args()
    return args

def main(args):
    ###########################
    ##                       ##
    ## Set local temp space  ##
    ##                       ##
    ###########################
    with tempfile.TemporaryDirectory() as tmpdirname:
        gs_loc = os.path.join(args.gs_bucket,os.path.basename(tmpdirname))
        task_fn= os.path.join(tmpdirname,'my-tasks.tsv')
        ###########################
        ##                       ##
        ## Write dsub task file  ##
        ##                       ##
        ###########################
        with open(task_fn,'w') as t_fh:
            print('--env CHUNK\t--input INFILE\t--output OUTFILE',file=t_fh)
            input_base = os.path.basename(args.input_data)
            input_link = os.path.join(gs_loc,input_base)
            output_base= os.path.basename(args.output_tag)
            output_link= os.path.join(gs_loc,output_base) + \
                         '__{}_' + '{}.bed'.format(args.job_count)
            for i in range(args.job_count):
                cmd_args = [str(i),input_link,output_link.format(i)]
                print('\t'.join(cmd_args),file=t_fh)
        #####################
        ##                 ##
        ## Send data to GS ##
        ##                 ##
        #####################
        gs_cmd = 'gsutil cp {} {}/'.format(args.input_data,gs_loc)
        gs_put = subprocess.run(gs_cmd.split())
        ###################################
        ##                               ##
        ## build dsub submission command ##
        ##                               ##
        ###################################
        vm_cmd = "python /app/hcr-ff/call_peaks.py ${INFILE} ${OUTFILE} " +\
                 "-ji ${CHUNK} " + "-jr {} ".format(args.job_count) +\
                 "-ws {} -ss {}".format(args.window_size, args.step_size)
        if args.no_offsets:
            vm_cmd += " --no_offsets"
        format_list = [args.billed_project, args.compute_zones, 
                       os.path.join(gs_loc,'logs'), task_fn, vm_cmd]
        dsub_cmd = "dsub --provider google-v2 --project {} --zones {} " +\
                   "--logging {} --machine-type n1-highmem-8 " +\
                   "--boot-disk-size 250 --disk-size 200 " +\
                   "--timeout 5h --tasks {} --image sjgosai/hff-kit:0.1.9 " +\
                   "--command '{}' --wait"
        if args.preemptable:
            dsub_cmd += " --preemptible --retries 5"
        else:
            dsub_cmd += " --retries 3"
        dsub_cmd = dsub_cmd.format( *format_list )
        print(dsub_cmd)
        ##############
        ##          ##
        ## Run dsub ##
        ##          ##
        ##############
        dsub_proc = subprocess.run(dsub_cmd, shell=True)
        ##################################
        ##                              ##
        ## Retrieve data from GS bucket ##
        ##                              ##
        ##################################
        gs_cmd = 'gsutil cp {} {}/'.format(output_link.format('*'),tmpdirname)
        gs_get = subprocess.run(gs_cmd.split())
        ##############
        ##          ##
        ## cat BEDs ##
        ##          ##
        ##############
        cat_cmd = ['cat'] + \
                  [ os.path.join( tmpdirname, 
                                    os.path.basename(output_link.format(i)) )
                    for i in range(args.job_count) ] + \
                  ['>',args.output_tag+'.bed']
        write_out = subprocess.run(" ".join(cat_cmd),shell=True)
        ########################
        ##                    ##
        ## Clean up GS bucket ##
        ##                    ##
        ########################
        gs_cmd = 'gsutil rm -r {}'.format(gs_loc)
        gs_rm  = subprocess.run(gs_cmd.split())

if __name__ == '__main__':
    args = get_args()
    main(args)
