/*
 * Run Trim Galore! on FastQ files
 */
process trimGalore {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }

    container 'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0'

    // Add a tag to identify the process
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_val_*.fq.gz") // Output trimmed FastQ files (R1 and R2 if paired-end)

    script:
    """
    echo "Running Trim Galore!"

    # Check the number of files in reads and run trim_galore accordingly
    if [ -f "${reads[0]}" ] && [ -f "${reads[1]}" ]; then
        trim_galore --paired ${reads[0]} ${reads[1]}
    elif [ -f "${reads[0]}" ]; then
        trim_galore ${reads[0]}
    else
        echo "No valid read files found for sample ${sample_id}"
        exit 1
    fi

    echo "Trim Galore! Complete"
    """
}