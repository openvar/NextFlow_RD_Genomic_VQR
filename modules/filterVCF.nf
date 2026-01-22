process filterVCF {
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$vcfFile"

    publishDir("$params.outdir/VCF", mode: "copy")

    input:
    tuple val(sample_id), file(vcfFile), file(vcfIndex)
    path indexFiles

    output:
    tuple val(sample_id), file("*_filtered.vcf")

    script:
    def isDegradedDNA = params.degraded_dna ? 'true' : 'false'
    """
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

    outputVcf="\$(basename ${vcfFile} .vcf)_filtered.vcf"

    if [ "$isDegradedDNA" == "true" ]; then
        echo "Running variant filtration for degraded DNA (2 x coverage)"
        gatk VariantFiltration -R "\${genomeFasta}" -V "${vcfFile}" -O "\${outputVcf}" \
            --filter-name "LowQUAL" --filter-expression "float(QUAL) < 150" \
            --filter-name "LowQD" --filter-expression "float(QD) < 2.0" \
            --filter-name "LowCoverage" --filter-expression "DP < 5" \
            --filter-name "HighFS" --filter-expression "float(FS) > 60.0" \
            --filter-name "HighSOR" --filter-expression "float(SOR) > 3.0" \
            --filter-name "LowMQ" --filter-expression "float(MQ) < 60.0" \
            --genotype-filter-name "LowGQ" --genotype-filter-expression "GQ < 30" \
            --set-filtered-genotype-to-no-call
    else
        echo "Running variant filtration for standard DNA (10x+ coverage)"
        gatk VariantFiltration -R "\${genomeFasta}" -V "${vcfFile}" -O "\${outputVcf}" \
            --filter-name "LowQUAL" --filter-expression "float(QUAL) < 500" \
            --filter-name "HighQD" --filter-expression "float(QD) < 2.0" \
            --filter-name "LowCoverage" --filter-expression "DP < 30" \
            --filter-name "HighFS" --filter-expression "float(FS) > 60.0" \
            --filter-name "HighSOR" --filter-expression "float(SOR) > 3.0" \
            --filter-name "LowMQ" --filter-expression "float(MQ) < 60.0" \
            --genotype-filter-name "LowGQ" --genotype-filter-expression "GQ < 30" \
            --set-filtered-genotype-to-no-call
    fi

    echo "Variant Filtering for Sample: ${vcfFile} Complete"
    """
}
