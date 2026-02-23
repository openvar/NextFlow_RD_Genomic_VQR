/*
 * Run fastp on the fastq files to trim adapters and low quality reads
 */

process fastp {

    container 'community.wave.seqera.io/library/fastp:1.1.0--08aa7c5662a30d57'

    // Add a tag to identify the process
    tag "$sample_id"

    // Specify the output directory for the fastp results
    publishDir("$params.outdir/fastp", mode: "copy")

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastp_${sample_id}_output/*"

    script:
    """
    echo "Running fastp"
    mkdir -p fastp_${sample_id}_output

    fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o fastp_${sample_id}_output/out.${sample_id}.R1.fastq.gz \
        -O fastp_${sample_id}_output/out.${sample_id}.R2.fastq.gz

    echo "fastp complete"
    """
}