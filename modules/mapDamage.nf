process mapDamage2 {

    label 'process_low'
    container 'quay.io/biocontainers/mapdamage2:2.2.2--pyr43hdfd78af_0'

    tag "$bamFile"

    // Publish mapDamage output to the specified directory
    publishDir("$params.outdir/mapDamage", mode: "copy")

    input:
    tuple val(sample_id), file(bamFile), file(baiFile)
    path indexFiles

    output:
    tuple val(sample_id), file("${bamFile.baseName}_rescaled.bam")

    script:
    """
    echo "Running mapDamage2 for sample ${sample_id}"

    # Define the genome file
    if [[ -n "${params.genome_file}" ]]; then
        genomeFasta=\$(basename "${params.genome_file}")
    else
        genomeFasta=\$(find -L . -name '*.fasta' | head -n 1)
    fi

    echo "Genome File: \${genomeFasta}"

    # Rename the dictionary file if it exists
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    # Define the output directory based on sample_id
    outputDir="${params.outdir}/mapDamage/results_${sample_id}_sorted"
    mkdir -p "\${outputDir}"

    # Extract base name for BAM file
    bamBaseName=\$(basename "${bamFile}" .bam)

    # Run mapDamage2 on the input BAM file, specifying output directory for rescaled BAM
    mapDamage -i "${bamFile}" \\
              -r "\${genomeFasta}" \\
              --rescale \\
              -d "\${outputDir}"

    # Move the rescaled BAM file to the work directory so it’s directly accessible by Nextflow
    mv "\${outputDir}/\${bamBaseName}.rescaled.bam" "./\${bamBaseName}_rescaled.bam"

    echo "mapDamage2 Complete for sample ${sample_id}"
    """
}
