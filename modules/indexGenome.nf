/*
 * Define the indexGenome process that creates a BWA index
 * given the genome fasta file
 */
process indexGenome {

    label 'process_low'
    container 'variantvalidator/indexgenome:1.1.0'


    // Publish indexed files to the specified directory
    publishDir("$params.outdir/GENOME_IDX", mode: "copy")

    input:
    path genomeFasta

    output:
    tuple path(genomeFasta), path("${genomeFasta}.*")

    script:
    """
    echo "Running Index Genome"

    # Generate BWA index
    bwa index "${genomeFasta}"

    # Generate samtools faidx
    samtools faidx "${genomeFasta}"

    # Generate Fasta dict
    picard CreateSequenceDictionary R="${genomeFasta}" O="${genomeFasta}.dict"

    echo "Genome Indexing complete."
    """
}
