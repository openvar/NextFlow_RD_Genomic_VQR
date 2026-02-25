/*
 * Run fastq on the read fastq files
 */
process FASTQC {

    label 'process_single'

    container 'variantvalidator/fastqc:0.12.1'

    // Add a tag to identify the process
    tag "${sample_id} (${stage})"

    // Specify the output directory for the FASTQC results
    publishDir("$params.outdir/FASTQC/${stage}", mode: "copy")

    input:
    tuple val(sample_id), val(stage), path(reads)

    output:
    tuple val(sample_id), val(stage), path("*.zip"), path("*.html")

    script:
    """
    echo "Running FASTQC"

    # Check the number of files in reads and run fastqc accordingly
    if [ -f "${reads[0]}" ] && [ -f "${reads[1]}" ]; then
        fastqc ${reads[0]} ${reads[1]}
    elif [ -f "${reads[0]}" ]; then
        fastqc ${reads[0]}
    else
        echo "No valid read files found for sample ${sample_id}"
        exit 1
    fi

    echo "FASTQC Complete"
    """
}
