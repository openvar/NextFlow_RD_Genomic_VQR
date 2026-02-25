/*
 * Align reads to the indexed genome using Minimap2
 */

process alignReadsMinimap2 {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }

    container 'quay.io/biocontainers/minimap2:2.28--h5bf99c6_0'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path genomeIndex

    output:
    tuple val(sample_id), file("${sample_id}.sam")

    script:
    """
    echo "Aligning reads with Minimap2 for sample ${sample_id}"

    read_args=${reads instanceof List ? reads.join(' ') : reads}

    # Align reads
    minimap2 -ax sr \
        # Add read group information
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
        $genomeIndex $read_args \
        > ${sample_id}.sam # Output SAM file

    echo "Minimap2 alignment complete for sample ${sample_id}"
    """
}

/*
 * Convert SAM to sorted BAM using samtools
 */

process samToSortedBam {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }

    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(samFile)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    echo "Converting SAM to sorted BAM for sample ${sample_id}"
    samtools sort -o ${sample_id}.bam ${samFile}
    rm ${samFile}
    echo "SAM to sorted BAM conversion complete for sample ${sample_id}"
    """
}