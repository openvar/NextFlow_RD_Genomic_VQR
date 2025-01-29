# Shell Commands

### Collapse reads

```shell
docker run --rm -v $(pwd):/input -v $(pwd):/output quay.io/biocontainers/adapterremoval:2.3.4--pl5321h6dccd9a_0 \
AdapterRemoval \
--file1 /input/sample_R1.fastq \
--file2 /input/sample_R2.fastq \
--basename /output/.fastq \
--qualitybase 33 \
--minlength 15 \
--trimqualities \
--minquality 15 \
--collapse
```

### Cat output files to include collapsed reads and unpaired reads to maximise data recovery
```shell
cat \
.fastq.singleton.truncated \
.fastq.collapsed \
> \
sample_combined.fastq
```

### Finally reomve any remaining adapter and trailing N bases
```shell
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--trim-n -q 20 -m 30 --cores=16 \
-o sample_trimmed_ca.fastq \
sample_combined.fastq
```

```
gunzip HF_combined_pretrimmed_ca.fastq && \
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--trim-n -q 20 -m 30 --cores=16 \
-o HF_combined_trimmed_ca.fastq \
HF_combined_pretrimmed_ca.fastq  && \
bgzip HF_combined_trimmed_ca.fastq

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--trim-n -q 20 -m 30 --cores=8 \
-o IMF_combined_trimmed_ca.fastq \
IMF_combined_pretrimmed_ca.fastq && \

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--trim-n -q 20 -m 30 --cores=8 \
-o IMF_HF_combined_trimmed_ca.fastq \
IMF_HF_combined_pretrimmed_ca.fastq && \
bgzip IMF_HF_combined_trimmed_ca.fastq
```