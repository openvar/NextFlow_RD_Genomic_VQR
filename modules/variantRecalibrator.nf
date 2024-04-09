process variantRecalibrator {

    label 'process_medium'
    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$vcf"

    // Publish VCF files to the specified directory
    publishDir("$params.outdir/VCF", mode: "copy")

    input:
    tuple val(sample_id), file(vcf)
    val knownSitesArgs
    path genome
    path qsrc_vcf

    output:
    tuple val(sample_id), file("${vcf.baseName}.recalibrated.vcf")

    script:
    def knownSitesArgsStr = knownSitesArgs.join(' ')
    """
    echo "Running VQSR"

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

    # Run VariantRecalibrator for SNPs
    gatk VariantRecalibrator \
        -R "\${genomeFasta}" \
        -V ${vcf} \
        ${knownSitesArgsStr} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
        -O ${vcf.baseName}.recalibrated_SNP.vcf

    # Apply VQSR for SNPs
    gatk ApplyVQSR \
        -R "\${genomeFasta}" \
        -V ${vcf} \
        --truth-sensitivity-filter-level 99.0 \
        -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
        -recal-file ${vcf.baseName}.recalibrated_SNP.recal \
        -mode SNP \
        -O ${vcf.baseName}.output_SNP.vcf

    # Run VariantRecalibrator for indels
    gatk VariantRecalibrator \
        -R "\${genomeFasta}" \
        -V ${vcf} \
        ${knownSitesArgsStr} \
        -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
        -mode INDEL \
        -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
        -O ${vcf.baseName}.recalibrated_INDEL.vcf

    # Apply VQSR for indels
    gatk ApplyVQSR \
        -R "\${genomeFasta}" \
        -V ${vcf} \
        --truth-sensitivity-filter-level 99.0 \
        -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
        -recal-file ${vcf.baseName}.recalibrated_INDEL.recal \
        -mode INDEL \
        -O ${vcf.baseName}.output_INDEL.vcf

    # Merge the recalibrated SNP and INDEL VCFs
    gatk MergeVcfs \
        -I ${vcf.baseName}.output_SNP.vcf \
        -I ${vcf.baseName}.output_INDEL.vcf \
        -O ${vcf.baseName}.recalibrated.vcf

    echo "VQSR Complete"
    """
}