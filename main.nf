// Use newest nextflow dsl
nextflow.enable.dsl = 2

// Print pipeline configuration
log.info """\
    ============================================
            Simple DNASeq Pipeline Configuration
    ============================================
    samplesheet     : ${params.samplesheet}
    genome          : ${params.genome_file}
    genome index    : ${params.genome_index_files}
    index genome    : ${params.index_genome}
    qsr truth vcfs  : ${params.qsrVcfs}
    output directory: ${params.outdir}
    fractions       : ${params.fractions}
""".stripIndent()

// Include modules in sequence
include { indexGenome } from './modules/indexGenome'
include { FASTQC } from './modules/FASTQC'
include { alignReads } from './modules/alignReads'
include { downsampleBam } from './modules/downsampleBam'
include { sortBam } from './modules/sortBam'
include { markDuplicates } from './modules/markDuplicates'
include { indexBam } from './modules/indexBam'
include { baseRecalibrator } from './modules/BQSR'
include { haplotypeCaller } from './modules/haplotypeCaller'
include { filterVCF } from './modules/filterVCF'

workflow {
    // Create a fractions channel
    fractions_ch = Channel.fromList(params.fractions)

    // User decides to index genome or not
    if (params.index_genome){
        // Flatten as is of format [fasta, [rest of files..]]
        indexed_genome_ch = indexGenome(params.genome_file).flatten()
        indexed_genome_ch.view()
    }
    else {
        indexed_genome_ch = Channel.fromPath(params.genome_index_files)
        indexed_genome_ch.view()
    }

    // Create qsrc_vcf_ch channel
    qsrc_vcf_ch = Channel.fromPath(params.qsrVcfs)
    qsrc_vcf_ch.view()

    // Set channel to gather read_pairs
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], [row[1], row[2]]) }
    read_pairs_ch.view()

    // Run FASTQC on read pairs
    FASTQC(read_pairs_ch)

    // Align reads to the indexed genome
    // Added index_genome_ch channel to make it wait for index if one not already created
    align_ch = alignReads(read_pairs_ch, indexed_genome_ch.collect())

    // Downsample BAM files
    // First combine fractions with bam channel
    // This means every downsample process will get its own instance so downsampling can run in parallel
    bam_fracs_ch = fractions_ch.combine(align_ch)
    downsample_ch = downsampleBam(bam_fracs_ch)

    // Sort BAM files
    // Removed collect to run sorting in parallel
    sort_ch = sortBam(downsample_ch)

    // Mark duplicates in BAM files
    // Removed collect to run sorting in parallel
    mark_ch = markDuplicates(sort_ch)

    // Index the BAM files and collect the output channel
    // Removed collect to run sorting in parallel
    // Changed indexBam to output sample_id, bam, bai file in tuple so that correct bai file stays with bamfile
    indexed_bam_ch = indexBam(mark_ch)

    // Create a channel from qsrVcfs. The process creates a text file with the paths to the vcf files
    knownSites_ch = Channel
        .fromPath(params.qsrVcfs)
        .filter { file -> file.getBaseName().endsWith('.vcf') }
        .map { file -> return "--known-sites " + file.getBaseName() }
        .collect()
    knownSites_ch.view()

    // Run BQSR on indexed BAM files
    bqsr_ch = baseRecalibrator(indexed_bam_ch, knownSites_ch, indexed_genome_ch.collect(), qsrc_vcf_ch.collect())

    // Call Variants
    // Removed collect from vcf bam ch
    // Altered haplotype caller to take in the index bam ch
    // Replaced params.genome_index_files with indexed_genome_ch
    vcf_call_ch = haplotypeCaller(bqsr_ch, indexed_genome_ch.collect())

    // Filter Variants
    // Removed collect from vcf channel
    // Replaced params.genome_index_files with indexed_genome_ch
    filterVCF(vcf_call_ch, indexed_genome_ch.collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nworkflow is done!\n" : "Oops .. something went wrong" )
}