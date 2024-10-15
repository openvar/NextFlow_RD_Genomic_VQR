// Use newest nextflow dsl
nextflow.enable.dsl = 2

process haplotypeCaller {
    label 'process_medium'
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

process mergeVCFs {
    label 'process_medium'
    container 'variantvalidator/gatk4:4.3.0.0'

    input:
    tuple val(sample_ids), path(gvcf_files), path(gvcf_index_files)

    output:
    tuple val(merged_sample_id), path("${merged_sample_id}.db")

    script:
    def merged_sample_id = sample_ids.join('_')
    def vcf_files_args = gvcf_files.collect { file -> "-V ${file}" }.join(' ')
    def chromosome_args = params.chromosomes_list.collect { "-L ${it}" }.join(' ')

    """
    echo "Merging VCFs for samples: ${gvcf_files.collect { it.baseName }.join(', ')}"

    gatk GenomicsDBImport --genomicsdb-workspace-path ${merged_sample_id}.db \
        ${vcf_files_args} \
        ${chromosome_args} \
        --reader-threads ${task.cpus}

    echo "VCF merge complete for samples: ${gvcf_files.collect { it.baseName }.join(', ')}"
    """
}

process GenotypeGVCFs {
    label 'process_medium'
    container 'variantvalidator/gatk4:4.3.0.0'

    input:
    tuple val(sample_id), path(gvcfFile) // Using path for merged GVCF files

    output:
    tuple val(sample_id), file("${sample_id}.vcf")

    script:
    def outputVcf = "${sample_id}.vcf"
    """
    echo "Genotyping GVCF for Sample: ${gvcfFile}"

    genomeFasta="\$(find -L . -name '*.fasta')"

    # Rename the dictionary file to the expected name if it exists
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    gatk GenotypeGVCFs -R "\${genomeFasta}" -V ${gvcfFile} -O ${outputVcf}

    echo "Genotyping complete: ${outputVcf}"
    """
}