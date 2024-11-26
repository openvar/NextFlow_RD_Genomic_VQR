process filterVCF {
    label 'process_low'
    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$vcfFile"

    // Publish VCF files to the specified directory
    publishDir("$params.outdir/VCF", mode: "copy")

    input:
    tuple val(sample_id), file(vcfFile)
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

    genomeFasta="\$(find -L . -name '*.fasta')"

    # Rename the dictionary file to the expected name
    mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"

    # Set output VCF filename with _filtered.vcf instead of .filtered.vcf
    outputVcf="\$(basename ${vcfFile} .vcf)_filtered.vcf"

    # If degraded DNA (3x coverage), use more relaxed filtering parameters, including MQ < 19 filter
    if [ "$isDegradedDNA" == "true" ]; then
        echo "Running variant filtration for degraded DNA (1x coverage)"
        gatk VariantFiltration -R "\${genomeFasta}" -V "${vcfFile}" -O "\${outputVcf}" \
            --filter-expression "QD < 2.0" --filter-name "LowQD" \
            --filter-expression "DP < 3" --filter-name "LowCoverage" \
            --filter-expression "FS > 60.0" --filter-name "HighFS" \
            --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
            --filter-expression "MQ < 19" --filter-name "LowMQ"

    # If standard DNA (10x or more coverage), use stricter parameters
    else
        echo "Running variant filtration for standard DNA (10x+ coverage)"
        gatk VariantFiltration -R "\${genomeFasta}" -V "${vcfFile}" -O "\${outputVcf}" \
            --filter-expression "QD < 2.0" --filter-name "LowQD" \
            --filter-expression "DP < 10" --filter-name "LowCoverage" \
            --filter-expression "FS > 60.0" --filter-name "HighFS" \
            --filter-expression "SOR > 3.0" --filter-name "HighSOR" \
            --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
            --filter-expression "MQRankSum < -12.5" --filter-name "LowMQRankSum" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "LowReadPosRankSum"
    fi

    # Print a message indicating the completion of variant filtration for the current sample
    echo "Variant Filtering for Sample: ${vcfFile} Complete"
    """
}
