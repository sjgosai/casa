# CASA
Analysis of Sabiti Lab HCR Flow-FISH data.

# Installation

`CASA` is not a Python library, but a collection of scripts to execute analysis. Dependencies are managed using `conda`, so if you don't have that start with:

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Once you have a recent version of `Anaconda` or `Miniconda`:

```
git clone https://github.com/sjgosai/casa.git
cd casa
conda env create -f casa_env.yml
```

Last, before using any code, activate your environment with:

```
conda activate casa
```

When everything is installed you should be ready to run all of the code in `./casa` and `./analysis`. However, `./casa/call_peaks.py` can process data in parallel on `GCP` using a docker environment with the above specs, and we've implemented a simple wrapper to do this which is dependent on `dsub`. To install:

# `GCP` and `dsub` setup

The easiest way to run `CASA` is using `GCP` and `dsub`. You can install `gsutil` and `dsub` anywhere (like on your MacBook or a VM) and run `CASA` on the cloud using `./src/wrap_peak_calling.py`. 

## Get a GCP account

If you don't have a `GCP` account, you can get a [free trial](https://cloud.google.com/free/) with $300 in credit with your gmail account. This should be more than enough to try out HCR analysis. Once you have an account, [create a billing project](https://console.cloud.google.com/projectcreate). Keep track of the project ID, you'll need it later.

## Install `gcloud` and `gsutil`

Next, setup the [Google Cloud SDK](https://cloud.google.com/deployment-manager/docs/step-by-step-guide/installation-and-setup "GCloud SDK Docs") and run `gcloud auth application-default login`. 

## Make a storage bucket for your data

Now, you need to configure where your input/output data will be stored in [Google Bucket Storage](https://cloud.google.com/storage/docs/quickstarts). For example, I want my input and output data to be stored in `gs://haddath/sgosai/hff/data/`. To do this, start by either using the `GCP` [console GUI](https://console.cloud.google.com/storage/browser) or `gsutil`:

```
gsutil mb -b on -l US gs://my-uniquely-named-bucket/
```

## Install `dsub`

Finally, you're ready to install dsub:

```
conda activate base
conda create --name dsub
pip install dsub
dsub --help
```


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
	--image sjgosai/casa-kit:0.2.1 \
	--command 'python /app/casa/call_peaks.py ${INFILE} ${OUTFILE} -ji ${CHUNK} -jr 20 -ws 100 -ss 100' \
	--wait &

```

Alternatively, we can use the batch job feature of `dsub` by specifying a task file:

my-tasks.tsv:
```
--env CHUNK	--input INFILE	--output OUTFILE
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
19	gs://haddath/sgosai/hff/data/FADS1_rep8detailed.txt	gs://haddath/sgosai/hff/data/FADS1_rep8__19_20.bed
```

And then following up with this command:

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
	--tasks my-tasks.tsv \
	--image sjgosai/casa-kit:0.2.1 \
	--command 'python /app/casa/call_peaks.py ${INFILE} ${OUTFILE} -ji ${CHUNK} -jr 20 -ws 100 -ss 100' \
	--wait &
```

Once this finishes running, you can pull the chunks from `bucket` storage and `cat` them together.
