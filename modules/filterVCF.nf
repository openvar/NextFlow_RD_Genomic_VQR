process filterVCF {
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$vcfFile"

    // Publish VCF files to the specified directory
    publishDir("$params.outdir/VCF", mode: "copy")

    input:
    tuple val(sample_id), file(vcfFile), file(vcfIndex)
    path indexFiles

    // Output channel for sample_id and filtered VCF files
    output:
    tuple val(sample_id), file("*_filtered.vcf")

    // Script section to run the process
    script:
    def isDegradedDNA = params.degraded_dna ? 'true' : 'false'
    """
    # Print a message indicating the start of the process for the current sample
    echo "Running Variant Filtration for Sample: ${vcfFile}"

    if [[ -n "${params.genome_file}" ]]; then
        genomeFasta=\$(basename ${params.genome_file})
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi

    echo "Genome File: \${genomeFasta}"

    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    # Set output VCF filename with _filtered.vcf instead of .filtered.vcf
    outputVcf="\$(basename ${vcfFile} .vcf)_filtered.vcf"

    # If degraded DNA (5x coverage), use more relaxed filtering parameters, including MQ < 40 filter
    if [ "$isDegradedDNA" == "true" ]; then
        echo "Running variant filtration for degraded DNA (2 x coverage)"
        gatk VariantFiltration -R "${genomeFasta}" -V "${vcfFile}" -O "${outputVcf}" \
            --filter-name "LowCoverage" --filter-expression "DP < 5" \
            --filter-name "HighFS" --filter-expression "FS > 60.0" \
            --filter-name "HighSOR" --filter-expression "SOR > 3.0" \
            --filter-name "LowMQ" --filter-expression "MQ < 40.0" \
            --genotype-filter-name "LowGQ" --genotype-filter-expression "GQ < 20" \
            --set-filtered-genotype-to-no-call \

    # If standard DNA use stricter parameters
    else
        echo "Running variant filtration for standard DNA (10x+ coverage)"
        gatk VariantFiltration -R "${genomeFasta}" -V "${vcfFile}" -O "${outputVcf}" \
            --filter-name "LowCoverage" --filter-expression "DP < 5" \
            --filter-name "HighFS" --filter-expression "FS > 60.0" \
            --filter-name "HighSOR" --filter-expression "SOR > 3.0" \
            --filter-name "LowMQ" --filter-expression "MQ < 60.0" \
            --filter-name "LowMQRankSum" --filter-expression "MQRankSum < -12.5" \
            --filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0"  \
            --genotype-filter-name "LowGQ" --genotype-filter-expression "GQ < 30" \
            --set-filtered-genotype-to-no-call
    fi


    # Print a message indicating the completion of variant filtration for the current sample
    echo "Variant Filtering for Sample: ${vcfFile} Complete"
    """
}
