source activate hff

cd /n/local/hcr-ff/analysis

python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS1_rep8detailed.txt FADS1_rep8__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 100 -gs

python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS1_rep9detailed.txt FADS1_rep9__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 100 -gs

python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS2_rep6detailed.txt FADS2_rep6__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 100 -gs

python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS2_rep8detailed.txt FADS2_rep8__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 100 -gs

python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS3_rep6detailed.txt FADS3_rep6__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 100 -gs

python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS3_rep8detailed.txt FADS3_rep8__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 100 -gs

