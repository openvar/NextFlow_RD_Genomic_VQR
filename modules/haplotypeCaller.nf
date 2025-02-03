// Use newest nextflow dsl
nextflow.enable.dsl = 2

process haplotypeCaller {
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }
    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$bamFile"

    input:
    tuple val(sample_id), file(bamFile), file(bamIndex)
    path indexFiles

    output:
    tuple val(sample_id), file("*.vcf"), file("*.vcf.idx")

    script:
    """
    echo "Running HaplotypeCaller for Sample: ${bamFile}"

    if [[ -n ${params.genome_file} ]]; then
        genomeFasta=\$(basename ${params.genome_file})
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi

    echo "Genome File: \${genomeFasta}"

    # Rename the dictionary file to the expected name if it exists
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    outputVcf="\$(basename ${bamFile} _sorted_dedup_recalibrated.bam).vcf"

    # Use GATK HaplotypeCaller to call variants in gVCF mode with specified annotations
    gatk HaplotypeCaller -R "\${genomeFasta}" -I ${bamFile} -O "\${outputVcf}" -ERC GVCF \
        -A BaseQuality -A DepthPerSampleHC -A MappingQuality -A QualByDepth \
        -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A StrandOddsRatio \
        -A MappingQualityZero -A InbreedingCoeff -A BaseQualityRankSumTest -A HaplotypeFilteringAnnotation

    echo "Sample: ${sample_id} VCF: \${outputVcf}"

    # Index the VCF file
    gatk IndexFeatureFile -I \${outputVcf}

    echo "Variant Calling for Sample: ${sample_id} Complete"
    """
}

