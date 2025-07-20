/*
 * Sort the BAM files and filter with samtools view
 */
process sortBam {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/indexgenome:1.1.0'

    tag "$bamFile"

    input:
    tuple val(sample_id), file(bamFile)

    output:
    tuple val(sample_id), file("${bamFile.baseName}_filtered_sorted.bam")

    script:
    """
    echo "Running Sort and Filter Bam"

    outputBam="\$(basename ${bamFile} .bam)_filtered_sorted.bam"

    # Use samtools to filter the BAM file with a minimum MAPQ score of 30, and then sort it
    samtools view -h -q 40 -b ${bamFile} | samtools sort -o \${outputBam} -

    echo "\${outputBam}"

    echo "BAM Sorting and Filtering complete"
    """
}
