process baseRecalibrator {

    label 'process_medium'
    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$bamFile"

    // Publish BQSR BAM files to the specified directory
    publishDir("$params.outdir/BAM", mode: "copy")

    input:
    tuple val(sample_id), file(bamFile), file(baiFile)
    val knownSites
    path indexFiles
    path qsrcVcfFiles

    output:
    tuple val(sample_id), file("${bamFile.baseName}_recalibrated.bam"), file("${bamFile.baseName}_recalibrated.bai")

    script:
    def knownSitesArgs = knownSites.join(' ')
    """
    echo "Running BQSR"

    if [[ -n params.genome_file ]]; then
        genomeFasta=\$(basename ${params.genome_file})
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi

    echo "Genome File: \${genomeFasta}"

    # Rename the dictionary file to the expected name if it exists
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    # Generate recalibration table for the input BAM file
    gatk --java-options "-Xmx8G" BaseRecalibrator \
        -R "\${genomeFasta}" \
        -I ${bamFile} \
        ${knownSitesArgs} \
        -O ${bamFile.baseName}.recal_data.table

    # Apply BQSR to the input BAM file
    gatk --java-options "-Xmx8G" ApplyBQSR \
        -R "\${genomeFasta}" \
        -I ${bamFile} \
        --bqsr-recal-file ${bamFile.baseName}.recal_data.table \
        -O ${bamFile.baseName}_recalibrated.bam

    # Index the recalibrated BAM file
    samtools index ${bamFile.baseName}_recalibrated.bam ${bamFile.baseName}_recalibrated.bai

    echo "BQSR Complete"
    """
}