# NextFlow_RD_Genomic

## Description

A simple base Rare disease and germline genomics pipeline to test the effects of down-sampling on variant calling

## Basic Overview
Using the NextFlow workflow software to run the following pipeline

## Setup
To run the pipeline, we need to obtain 

- A genome build (GRCh38) - provided by the Broad institute - ***REQUIRED***
```bash
$ cd data/genome
$ wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
```
- FastQ sample (for workflow development) - ***REQUIRED***
```bash
$ cd ../samples
$ wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518253/SRR1518253_1.fastq.gz && \
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR151/003/SRR1518253/SRR1518253_2.fastq.gz
```

- FASTQ for low coverage WGS - ***OPTIONAL***
```
# Download the files using HTTPS
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR063/SRR063269/SRR063269_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR063/SRR063279/SRR063279_2.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR063/SRR063269/SRR063269_2.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR063/SRR063277/SRR063277_2.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR063/SRR063277/SRR063277_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR063/SRR063279/SRR063279_1.fastq.gz

# Concatenate the R1 files into NA19750_R1
cat SRR063269_1.fastq.gz SRR063277_1.fastq.gz SRR063279_1.fastq.gz > NA19750_lc_wgs_R1.fastq.gz

# Concatenate the R2 files into NA19750_R2
cat SRR063269_2.fastq.gz SRR063277_2.fastq.gz SRR063279_2.fastq.gz > NA19750_lc_wgs_R2.fastq.gz
```

- FASTQ for NA19750 Exome - ***OPTIONAL***
```
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR071/SRR071178/SRR071178_1.fastq.gz
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR071/SRR071178/SRR071178_2.fastq.gz
```

## For VQSR and BQSR download the following files and indexes into the relevant directory - ***REQUIRED***
```bash
# 1000G Omni 2.5
wget -L "https://zenodo.org/records/18267841/files/1000G_omni2.5.hg38.vcf.gz?download=1"
wget -L "https://zenodo.org/records/18267841/files/1000G_omni2.5.hg38.vcf.gz.tbi?download=1"

# 1000G Phase 1 high-confidence SNPs
wget -L "https://zenodo.org/records/18267841/files/1000G_phase1.snps.high_confidence.hg38.vcf.gz?download=1"
wget -L "https://zenodo.org/records/18267841/files/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi?download=1"

# HapMap 3.3
wget -L "https://zenodo.org/records/18267841/files/hapmap_3.3.hg38.vcf.gz?download=1"
wget -L "https://zenodo.org/records/18267841/files/hapmap_3.3.hg38.vcf.gz.tbi?download=1"

# dbSNP 138 (uncompressed VCF)
wget -L "https://zenodo.org/records/18267841/files/Homo_sapiens_assembly38.dbsnp138.vcf?download=1"
wget -L "https://zenodo.org/records/18267841/files/Homo_sapiens_assembly38.dbsnp138.vcf.idx?download=1"

# Mills + 1000G gold standard indels
wget -L "https://zenodo.org/records/18267841/files/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz?download=1"
wget -L "https://zenodo.org/records/18267841/files/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi?download=1"
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


## GRCh37 Download
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/p13.plusMT/hg19.p13.plusMT.fa.gz -O Homo_sapiens_assembly37.fasta.gz
```

## Syphilis genome
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/745/GCA_000008745.1_ASM874v1/GCA_000008745.1_ASM874v1_genomic.fna.gz -O syphilis.fa.gz
```

## Concatenate gzipped files
```bash
zcat file1.gz file2.gz | gzip > concatenated_file.gz
```
