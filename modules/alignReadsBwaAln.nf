process alignReadsBwaAln {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_long'
    }
    container 'variantvalidator/indexgenome:1.1.0'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)   // reads is a tuple of paths for paired-end reads
    path requiredIndexFiles

    output:
    tuple val(sample_id), file("${sample_id}.bam")

    script:
    """
    # Find the BWA index
    INDEX=\$(find -L ./ -name "*.amb" | sed 's/\\.amb\$//')

    echo "Running Align Reads with BWA-ALN"
    echo "\$INDEX"

    # Check if the input FASTQ files exist
    if [ -f "${reads[0]}" ]; then
        if [ -f "${reads[1]}" ]; then
            # Paired-end mode
            bwa aln -t ${task.cpus} \$INDEX ${reads[0]} > ${sample_id}_1.sai
            bwa aln -t ${task.cpus} \$INDEX ${reads[1]} > ${sample_id}_2.sai

            # Use BWA-SAMPE to align paired-end reads and pipe to samtools to generate BAM file
            bwa sampe \$INDEX ${sample_id}_1.sai ${sample_id}_2.sai ${reads[0]} ${reads[1]} |
            samtools view -b - |
            # Set MAPQ to 0 for unmapped reads using samtools
            awk '\$5 == 0 { \$5 = 0 } { print }' |
            samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam
        else
            # Single FASTQ mode
            bwa aln -t ${task.cpus} \$INDEX ${reads[0]} > ${sample_id}_1.sai
            bwa samse \$INDEX ${sample_id}_1.sai ${reads[0]} |
            samtools view -b - |
            # Set MAPQ to 0 for unmapped reads using samtools
            awk '\$5 == 0 { \$5 = 0 } { print }' |
            samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam
        fi
    else
        echo "Error: Read file ${reads[0]} does not exist for sample ${sample_id}."
        exit 1
    fi

    echo "Alignment complete"
    """
}
