/*
 * Align reads to the indexed genome
 */
process alignReadsBwaAln {

    label 'process_medium'
    container 'variantvalidator/indexgenome:1.1.0'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path requiredIndexFiles

    output:
    tuple val(sample_id), file("${sample_id}.bam")

    script:
    """
    # Find the BWA index
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    echo "Running Align Reads with BWA-ALN"
    echo "\$INDEX"

    # Step 1: Run BWA-ALN to generate .sai files
    bwa aln -t ${task.cpus} \$INDEX ${reads[0]} > ${sample_id}_1.sai
    bwa aln -t ${task.cpus} \$INDEX ${reads[1]} > ${sample_id}_2.sai

    # Step 2: Use BWA-SAMPE to align paired-end reads and pipe to samtools to generate BAM file
    bwa sampe \$INDEX ${sample_id}_1.sai ${sample_id}_2.sai ${reads[0]} ${reads[1]} |
    samtools view -b - |
    samtools addreplacerg -r "@RG\\tID:${sample_id}_1.fastq\\tSM:${sample_id}_2.fastq\\tPL:illumina" - |
    samtools addreplacerg -r "@RG\\tID:${sample_id}_2.fastq\\tSM:${sample_id}_1.fastq\\tPL:illumina" - > ${sample_id}.bam

    echo "Alignment complete"
    """
}
