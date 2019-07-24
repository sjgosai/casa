# hcr-ff
Analysis of Sabiti Lab HCR Flow-FISH data.

# `dsub` pipeline
To start, make sure you've setup `gcloud`, the [Google Cloud SDK](https://cloud.google.com/deployment-manager/docs/step-by-step-guide/installation-and-setup "GCloud SDK Docs") and run `gcloud auth application-default login`. Once you've done this, you need decide where your input/output data will be stored in Google Bucket Storage. For example, I want my input and output data to be stored in `gs://haddath/sgosai/hff/data/`.

## Submitting jobs
The peak caller uses 8 threads for computation and is executed on chunks of data. The user specifies which chunk each instance of the script runs on, and the total number of chunks. We can submit a job to process a chunk as follows:
```
dsub \
	--provider google-v2 \
	--project sabeti-encode \
	--zones "us-*" \
	--logging gs://haddath/sgosai/hff/logs \
	--machine-type n1-standard-8 \
	--boot-disk-size 250 \
	--disk-size 50 \
	--preemptible \
	--retries 3 \
	--env CHUNK=0
	--input INFILE=gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt \
	--output OUTFILE=gs://haddath/sgosai/hff/data/FADS1_rep8__0_20.bed \
	--image sjgosai/hff-kit:0.1.1 \
	--command 'python /app/hcr-ff/call_peaks.py ${INFILE} ${OUTFILE} -ji ${CHUNK} -jr 20 -ws 100 -ss 100' \
	--wait &

```

Alternatively, we can use the batch job feature of `dsub` by specifying a task file:

my-tasks.tsv
```--env CHUNK	--input INFILE	--output OUTFILE
0	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__0_20.bed
1	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__1_20.bed
2	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__2_20.bed
3	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__3_20.bed
4	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__4_20.bed
5	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__5_20.bed
6	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__6_20.bed
7	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__7_20.bed
8	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__8_20.bed
9	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__9_20.bed
10	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__10_20.bed
11	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__11_20.bed
12	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__12_20.bed
13	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__13_20.bed
14	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__14_20.bed
15	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__15_20.bed
16	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__16_20.bed
17	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__17_20.bed
18	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__18_20.bed
19	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__19_20.bed```

```dsub \
	--provider google-v2 \
	--project sabeti-encode \
	--zones "us-*" \
	--logging gs://haddath/sgosai/hff/logs \
	--machine-type n1-standard-8 \
	--boot-disk-size 250 \
	--disk-size 50 \
	--preemptible \
	--retries 3 \
	--tasks my-tasks.tsv \
	--image sjgosai/hff-kit:0.1.2 \
	--command 'python /app/hcr-ff/call_peaks.py ${INFILE} ${OUTFILE} -ji ${CHUNK} -jr 20 -ws 100 -ss 100' \
	--wait &```
