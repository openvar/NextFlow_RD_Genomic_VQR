/*
 * Align reads to the indexed genome
 */
process alignReadsBwaMem {

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
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    
    echo "Running Align Reads"

    echo "\$INDEX"

    # Align reads using BWA and generate BAM file
    bwa mem -M -B 3 -O 6 -E 1 -k 17 -w 120 -T 20 -t ${task.cpus} \$INDEX ${reads[0]} ${reads[1]} |
    samtools view -b - |
    samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam
    echo "Alignment complete"
    """
}
