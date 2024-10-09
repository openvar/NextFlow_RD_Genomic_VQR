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
include { alignReadsBwaMem } from './modules/alignReadsBwaMem'
include { downsampleBam } from './modules/downsampleBam'
include { sortBam } from './modules/sortBam'
include { markDuplicates } from './modules/markDuplicates'
include { indexBam } from './modules/indexBam'
include { baseRecalibrator } from './modules/BQSR'
include { haplotypeCaller } from './modules/haplotypeCaller'
include { variantRecalibrator } from './modules/variantRecalibrator'


workflow {
    // Create a fractions channel
    fractions_ch = Channel.fromList(params.fractions)

    // User decides to index genome or not
    if (params.index_genome){
        // Flatten as is of format [fasta, [rest of files..]]
        indexed_genome_ch = indexGenome(params.genome_file).flatten()
    }
    else {
        indexed_genome_ch = Channel.fromPath(params.genome_index_files)
    }

    // Create qsrc_vcf_ch channel
    qsrc_vcf_ch = Channel.fromPath(params.qsrVcfs)

    // Set channel to gather read_pairs
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], [row[1], row[2]]) }
    read_pairs_ch.view()

    // Run FASTQC on read pairs
    FASTQC(read_pairs_ch)

    // Align reads to the indexed genome
    align_ch = alignReadsBwaMem(read_pairs_ch, indexed_genome_ch.collect())
    align_ch.view()

    // Conditionally downsample BAM files
    if (params.downsample_bam) {
        // Combine fractions with bam channel
        bam_fracs_ch = fractions_ch.combine(align_ch)
        downsample_ch = downsampleBam(bam_fracs_ch)

        // Sort BAM files
        sort_ch = sortBam(downsample_ch)
    } else {
        // If downsampling is disabled, directly sort the aligned BAM files
        sort_ch = align_ch.map { bamFile ->
            def sample_id = bamFile.getBaseName().split('_')[0]
            return tuple(sample_id, bamFile)
        }
        sort_ch = sortBam(sort_ch)
    }

    // Mark duplicates in BAM files
    sort_ch.view()
    mark_ch = markDuplicates(sort_ch)

    // Index the BAM files and collect the output channel
    indexed_bam_ch = indexBam(mark_ch)

    // Create a channel from qsrVcfs
    knownSites_ch = Channel
        .fromPath(params.qsrVcfs)
        .filter { file -> file.getName().endsWith('.vcf.idx') }
        .map { file -> return "--known-sites " + file.getBaseName() }
        .collect()
    knownSites_ch.view()

    // Run BQSR on indexed BAM files
    bqsr_ch = baseRecalibrator(indexed_bam_ch, knownSites_ch, indexed_genome_ch.collect(), qsrc_vcf_ch.collect())

    // Call Variants
    vcf_call_ch = haplotypeCaller(bqsr_ch, indexed_genome_ch.collect())

    // Define a map of VCF files to resource options
    def resourceOptions = [
        'Homo_sapiens_assembly38.known_indels': '--resource:dbsnp_indels,known=true,training=false,truth=false,prior=15.0',
        'hapmap_3.3.hg38': '--resource:hapmap,known=false,training=false,truth=true,prior=15.0',
        '1000G_omni2.5.hg38': '--resource:omni,known=false,training=true,truth=false,prior=12.0',
        '1000G_phase1.snps.high_confidence.hg38': '--resource:1000G,known=true,training=true,truth=true,prior=10.0',
        'Homo_sapiens_assembly38.dbsnp138': '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0'
    ]

    knownSitesArgs_ch = Channel
        .fromPath(params.qsrVcfs)
        .map { file -> return file }
        .filter { file -> file.getName().endsWith('.vcf') }
        .map { file ->
            def basename = file.getBaseName()
            if (resourceOptions.containsKey(basename)) {
                return resourceOptions[basename.replace('.vcf', '')] + " ./" + basename + ".vcf"
            } else {
                return null
            }
        }
        .filter { it != null }
        .collect()
    knownSitesArgs_ch.view()

    // Apply Variant Recalibration
    variantRecalibrator(vcf_call_ch, knownSitesArgs_ch, indexed_genome_ch.collect(), qsrc_vcf_ch.collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nworkflow is done!\n" : "Oops .. something went wrong" )
}