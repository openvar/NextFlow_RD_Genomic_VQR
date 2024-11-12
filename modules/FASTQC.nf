/*
 * Run fastq on the read fastq files
 */
process FASTQC {

    label 'process_single'

    container 'variantvalidator/fastqc:0.12.1'

    // Add a tag to identify the process
    tag "$sample_id"

    // Specify the output directory for the FASTQC results
    publishDir("$params.outdir/FASTQC", mode: "copy")

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs/*"

    script:
    """
    echo "Running FASTQC"
    mkdir -p fastqc_${sample_id}_logs

    # Check the number of files in reads and run fastqc accordingly
    if [ -f "${reads[0]}" ] && [ -f "${reads[1]}" ]; then
        fastqc ${reads[0]} ${reads[1]} -o fastqc_${sample_id}_logs
    elif [ -f "${reads[0]}" ]; then
        fastqc ${reads[0]} -o fastqc_${sample_id}_logs
    else
        echo "No valid read files found for sample ${sample_id}"
        exit 1
    fi

    echo "FASTQC Complete"
    """
}
