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
    aligner         : ${params.aligner}
    read_trim       : ${params.read_trimming}
    variant_recalibration: ${params.variant_recalibration}
    ============================================
""".stripIndent()

// Include modules in sequence
include { indexGenome } from './modules/indexGenome'
include { leehom } from './modules/leehom'
include { FASTQC } from './modules/FASTQC'
include { sortBam } from './modules/sortBam'
include { markDuplicates } from './modules/markDuplicates'
include { indexBam } from './modules/indexBam'
include { baseRecalibrator } from './modules/BQSR'
include { haplotypeCaller } from './modules/haplotypeCaller'
include { mergeVCFs } from './modules/haplotypeCaller'
include { GenotypeGVCFs } from './modules/haplotypeCaller' // Include GenotypeGVCFs process
include { variantRecalibrator } from './modules/variantRecalibrator'

// Conditionally include the alignment process based on the aligner parameter
if (params.aligner == 'bwa-mem') {
    include { alignReadsBwaMem } from './modules/alignReadsBwaMem'
} else if (params.aligner == 'bwa-aln') {
    include { alignReadsBwaAln } from './modules/alignReadsBwaAln'
} else {
    error "Unsupported aligner: ${params.aligner}. Please specify 'bwa-mem' or 'bwa-aln'."
}

workflow {
    // Create channels and run processes in sequence
    indexed_genome_ch = params.index_genome ? indexGenome(params.genome_file).flatten() : Channel.fromPath(params.genome_index_files)
    qsrc_vcf_ch = Channel.fromPath(params.qsrVcfs)
    read_pairs_ch = Channel.fromPath(params.samplesheet).splitCsv(sep: '\t').map { row -> tuple(row[0], [row[1], row[2]]) }
    read_pairs_ch.view()

    // Conditionally run leehom if read trimming is enabled
    if (params.read_trimming) {
        trimmed_reads_ch = leehom(read_pairs_ch)
        trimmed_reads_ch.view()
    } else {
        trimmed_reads_ch = read_pairs_ch
    }

    // Run FASTQC on read pairs
    FASTQC(trimmed_reads_ch)

    // Align reads to the indexed genome
    if (params.aligner == 'bwa-mem') {
        align_ch = alignReadsBwaMem(trimmed_reads_ch, indexed_genome_ch.collect())
    } else if (params.aligner == 'bwa-aln') {
        align_ch = alignReadsBwaAln(trimmed_reads_ch, indexed_genome_ch.collect())
    }

    // Sort BAM files
    sort_ch = sortBam(align_ch)

    // Mark duplicates in BAM files
    mark_ch = markDuplicates(sort_ch)

    // Index the BAM files and collect the output channel
    indexed_bam_ch = indexBam(mark_ch)

    // Create a channel from qsrVcfs
    knownSites_ch = Channel.fromPath(params.qsrVcfs).filter { file -> file.getName().endsWith('.vcf.idx') }.map { file -> "--known-sites " + file.getBaseName() }.collect()

    // Run BQSR on indexed BAM files
    bqsr_ch = baseRecalibrator(indexed_bam_ch, knownSites_ch, indexed_genome_ch.collect(), qsrc_vcf_ch.collect())

    // Call Variants in gVCF mode and collect the outputs
    gvcf_ch = haplotypeCaller(bqsr_ch, indexed_genome_ch.collect())
        .collect { tuple ->
            [
                tuple[0],    // sample_id (first element)
                tuple[1],    // vcf_file (second element)
                tuple[2]     // vcf_index_file (third element)
            ]
        }
        .view()

    // Now we map to create separate lists for sample IDs, VCF files, and index files
    all_gvcf_ch = gvcf_ch
        .collect { listOfTuples ->
            def sample_ids = listOfTuples.collect { it[0] }  // Collect sample IDs
            def vcf_files = listOfTuples.collect { it[1] }   // Collect VCF files
            def vcf_index_files = listOfTuples.collect { it[2] } // Collect VCF index files
            return tuple(sample_ids, vcf_files, vcf_index_files)
        }

    // Merge VCFs using GenomicsDBImport
    merged_vcf_db = mergeVCFs(all_gvcf_ch)

    // Genotype merged GVCFs
    final_vcf_ch = merged_vcf_db.flatMap { db_tuple ->
        def merged_sample_id = db_tuple[0] // Merged sample ID
        def merged_db_path = db_tuple[1] // Merged DB path
        return GenotypeGVCFs(merged_sample_id, merged_db_path) // Genotype the merged GVCF
    }
    final_vcf_ch.view() // View final VCF files

    // Conditionally apply variant recalibration or filtering
    if (params.variant_recalibration) {
        // Define a map of VCF files to resource options
        def resourceOptions = [
            'Homo_sapiens_assembly38.known_indels': '--resource:dbsnp_indels,known=true,training=false,truth=false,prior=15.0',
            'hapmap_3.3.hg38': '--resource:hapmap,known=false,training=false,truth=true,prior=15.0',
            '1000G_omni2.5.hg38': '--resource:omni,known=false,training=true,truth=false,prior=12.0',
            '1000G_phase1.snps.high_confidence.hg38': '--resource:1000G,known=true,training=true,truth=true,prior=10.0',
            'Homo_sapiens_assembly38.dbsnp138': '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0'
        ]

        knownSitesArgs_ch = Channel.fromPath(params.qsrVcfs).map { file -> file }.filter { file -> file.getName().endsWith('.vcf') }.map { file -> resourceOptions[file.getBaseName().replace('.vcf', '')] + " ./" + file.getBaseName() + ".vcf" }.filter { it != null }.collect()
        knownSitesArgs_ch.view()

        filtered_vcf_ch = variantRecalibrator(final_vcf_ch, knownSitesArgs_ch, indexed_genome_ch.collect(), qsrc_vcf_ch.collect())
    } else {
        filtered_vcf_ch = filterVCF(final_vcf_ch)
    }
    filtered_vcf_ch.view() // View final VCF files after filtering or recalibration
}

workflow FASTQC_only {
    // Create channels and run processes in sequence
    read_pairs_ch = Channel.fromPath(params.samplesheet).splitCsv(sep: '\t').map { row -> tuple(row[0], [row[1], row[2]]) }
    read_pairs_ch.view()

    // Run FASTQC on read pairs
    FASTQC(read_pairs_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nworkflow is done!\n" : "Oops .. something went wrong" )
}

