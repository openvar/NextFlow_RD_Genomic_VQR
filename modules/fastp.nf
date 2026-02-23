/*
 * Run fastp on the fastq files to trim adapters and low quality reads
 */

process fastp {

    container 'community.wave.seqera.io/library/fastp:1.1.0--08aa7c5662a30d57'

    input:
    path input_fastqs

    output:
}