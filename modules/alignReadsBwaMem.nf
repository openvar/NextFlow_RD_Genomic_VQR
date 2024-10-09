/*
 * Align reads to the indexed genome
 */
process alignReadsBwaMem {

    label 'process_medium'
    container 'variantvalidator/indexgenome:1.1.0'
     
    tag "alignReadsBwaMem"
    
    input:
    tuple val(sample_id), path(reads)
    path requiredIndexFiles

    output:
    file("${sample_id}.bam")

    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    
    echo "Running Align Reads"

    echo "\$INDEX"

    # Align reads using BWA and generate BAM file
    bwa mem -t ${task.cpus} \$INDEX ${reads[0]} ${reads[1]} |
    samtools view -b - |
    samtools addreplacerg -r "@RG\\tID:${sample_id}_1.fastq\\tSM:${sample_id}_2.fastq\\tPL:illumina" - |
    samtools addreplacerg -r "@RG\\tID:${sample_id}_2.fastq\\tSM:${sample_id}_1.fastq\\tPL:illumina" - > ${sample_id}.bam
    echo "Alignment complete"
    """
}
