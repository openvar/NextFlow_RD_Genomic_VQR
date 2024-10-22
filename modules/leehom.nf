/*
 * Run leehom on the input FASTQ files
 */
process leehom {

    label 'process_medium'
    container 'variantvalidator/leehom:0.0.1'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), file("${sample_id}_leehom.fastq")

    script:
    """
    echo "Running leehom for Sample: ${sample_id}"

    outputFastq="${sample_id}_leehom.fastq"

    # Check if degraded DNA flag is set and adjust leehom command accordingly
    if [ "${params.degraded_dna}" == "true" ]; then
        echo "Using aDNA enhancements for degraded DNA."
        leehom --no-align --ancient --fragadapter -i ${reads[0]} -I ${reads[1]} -o \${outputFastq}
    else
        echo "Running leehom without aDNA enhancements."
        leehom --no-align -i ${reads[0]} -I ${reads[1]} -o \${outputFastq}
    fi

    echo "Sample: ${sample_id} leehom output: \${outputFastq}"

    echo "leehom processing for Sample: ${sample_id} Complete"
    """
}
