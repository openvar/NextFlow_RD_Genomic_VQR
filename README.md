# NextFlow_RD_Genomic

# Tasks

- Create the directory data/qsr_vcfs
- Populate with https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
- Populate with https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
- Make the modules/BSQR.nf module
- Search main.nf for changes (note the channel creation and usage)
- Update the nextflow.config

## Description

A simple base Rare disease and germline genomics pipeline to test the effects of down-sampling on variant calling

## Basic Overview
Using the NextFlow workflow software to run the following pipeline

### Pipeline
Index genome > Fastqc analysis > Align reads > Downsample bam files > Sort bam > Mark duplicates > Index bam > 
Call variants > Hard filter

## Setup
To run the pipeline, we need to obtain 

- A genome build (GRCh38) - provided by the Broad institute
```bash
$ cd data/genome
$ wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
```
- FastQ sample (for workflow development)
```bash
$ cd ../samples
$ wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518253/SRR1518253_1.fastq.gz && \
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518253/SRR1518253_2.fastq.gz && \
```
- When scaling up, FastQ samples
```bash
$ wget https://genomics.viapath.co.uk/benchmark/files/FASTQ/NA12878_WES.zip
```
Note: the file will need to be unzipped and compressed using bgzip

## For VQSR and BQSR download the following files and indexes into the relevant directory
```bash
$ wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz &&
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz &&
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz &&
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz &&
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx &&
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi &&
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi &&
wget https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
```

## Samplesheet
For paired end alignment use the format

sampleID	fastq_path	Sex
```
sample1 /path/to/sample1_1.fastq.gz	/path/to/sample1_2.fastq.gz	1
sample2 /path/to/sample2_1.fastq.gz	/path/to/sample2_2.fastq.gz	2
```

and for single end or collapsed reads use

sampleID	fastq_path	Sex
```
sample1 /path/to/sample1.fastq.gz	1
sample2 /path/to/sample2.fastq.gz	2
```

## Running the pipeline
```bash
# Using Docker
$ nextflow run -profile docker main.nf

# Using docker in singularity
$ nextflow run -profile singularity main.nf
```
Note: Refer to the nexrflow.config and nextflow_schema.json for parameter selection. 

## Validating the pipeline
See [https://genomics.viapath.co.uk/benchmark](https://genomics.viapath.co.uk/benchmark)

## DNANexus applet setup (A local applet for basic testing)
- DNANexus Python Bindings [Documentation](https://github.com/dnanexus/dx-toolkit) 
- [Install the app](https://documentation.dnanexus.com/downloads) 
```bash
pip install -r requirements.txt
```
- Routine maintenance
Periodically update dxpy
```bash
$ pip install --upgrade dxpy
```

### DNANexus Tutorial
- [Overview videos](https://documentation.dnanexus.com/getting-started)
- [Developer tutorial](https://documentation.dnanexus.com/getting-started/developer-quickstart)
```bash
$ dx select <your-project-name>
$ dx build --nextflow
```