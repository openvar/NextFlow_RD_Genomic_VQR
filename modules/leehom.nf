/*
 * Run leehom on the input FASTQ files
 */
process leehom {

    label 'process_medium'
    container 'quay.io/biocontainers/leehom:1.2.15--hdc46a4b_6'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), file("${sample_id}_leehom.fastq")

    script:
    """
    echo "Running leehom for Sample: ${sample_id}"

    outputFastq="${sample_id}_leehom.fastq"

    # Run leehom on the input FASTQ files
    leehom --no-align -i ${reads[0]} -I ${reads[1]} -o \${outputFastq}

    echo "Sample: ${sample_id} leehom output: \${outputFastq}"

    echo "leehom processing for Sample: ${sample_id} Complete"
    """
}