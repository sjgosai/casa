source activate casa

cd /n/local/casa/analysis

#python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS1_rep8detailed.txt FADS1_rep8__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 50 -gs

#python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS1_rep9detailed.txt FADS1_rep9__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 50 -gs

#python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS2_rep8detailed.txt FADS2_rep8__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 50 -gs

python /n/local/casa/casa/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS2_rep11detailed.txt FADS2_rep11__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 50 -gs

#python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS3_rep7detailed.txt FADS3_rep7__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 50 -gs

#python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FADS3_rep8detailed.txt FADS3_rep8__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 50 -gs

#python /n/local/hcr-ff/hcr-ff/call_peaks.py /n/local/hcr-ff/data/FASTQ/FEN1_rep11detailed.txt FEN1_rep11__$1_20.bed -ji $1 -jr 20 -ws 100 -ss 50 -gs
