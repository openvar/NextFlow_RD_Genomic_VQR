/*
 * Sort the BAM files
 */
process sortBam {

    label 'process_low'
    container 'variantvalidator/indexgenome:1.1.0'

    tag "$bamFile"

    input:
    tuple val(sample_id), file(bamFile)

    output:
    tuple val(sample_id), file("${bamFile.baseName}_sorted.bam")

    script:
    """
    echo "Running Sort Bam"

    outputBam="\$(basename ${bamFile} .bam)_sorted.bam"

    # Use samtools to sort the input BAM file and save the output to the sorted BAM file
    samtools sort -o \${outputBam} ${bamFile}

    echo "\${outputBam}"

    echo "BAM Sorting complete"
    """
}
