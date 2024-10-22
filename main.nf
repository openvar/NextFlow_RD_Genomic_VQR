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
    degraded_dna    : ${params.degraded_dna}
    variant_recalibration: ${params.variant_recalibration}
    identity_analysis: ${params.identity_analysis}
    ============================================
""".stripIndent()

// Conditionally include modules
if (params.index_genome) {
    include { indexGenome } from './modules/indexGenome'
}
if (params.read_trimming) {
    include { leehom } from './modules/leehom'
}
include { FASTQC } from './modules/FASTQC'
include { sortBam } from './modules/sortBam'
include { markDuplicates } from './modules/markDuplicates'
include { indexBam } from './modules/indexBam'
include { baseRecalibrator } from './modules/BQSR'
include { haplotypeCaller } from './modules/haplotypeCaller'
include { combineGVCFs } from './modules/haplotypeCaller'
include { genotypeGVCFs } from './modules/haplotypeCaller'
if (params.variant_recalibration) {
    include { variantRecalibrator } from './modules/variantRecalibrator'
} else {
    include { filterVCF } from './modules/filterVCF'
}
if (params.identity_analysis) {
    include { identityAnalysis } from './modules/identityAnalysis'
}

// Conditionally include the alignment process based on the aligner parameter
if (params.aligner == 'bwa-mem') {
    include { alignReadsBwaMem } from './modules/alignReadsBwaMem'
} else if (params.aligner == 'bwa-aln') {
    include { alignReadsBwaAln } from './modules/alignReadsBwaAln'
} else {
    error "Unsupported aligner: ${params.aligner}. Please specify 'bwa-mem' or 'bwa-aln'."
}

workflow {
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

    // Run HaplotypeCaller on BQSR files
    gvcf_ch = haplotypeCaller(bqsr_ch, indexed_genome_ch.collect())
        .collect()

    // Now we map to create separate lists for sample IDs, VCF files, and index files
    all_gvcf_ch = gvcf_ch
        .collect { listOfTuples ->
            def sample_ids = listOfTuples.collate(3).collect { it[0] }   // Collect sample IDs from every 3rd element
            def vcf_files = listOfTuples.collate(3).collect { it[1] }    // Collect VCF files
            def vcf_index_files = listOfTuples.collate(3).collect { it[2] } // Collect VCF index files
            return tuple(sample_ids, vcf_files, vcf_index_files)
        }

    // Combine GVCFs
    combined_gvcf_ch = combineGVCFs(all_gvcf_ch, indexed_genome_ch.collect())

    // Run GenotypeGVCFs
    final_vcf_ch = genotypeGVCFs(combined_gvcf_ch, indexed_genome_ch.collect())

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

        knownSitesArgs_ch = Channel.fromPath(params.qsrVcfs).map { file -> file }.filter { file.getName().endsWith('.vcf') }.map { file -> resourceOptions[file.getBaseName().replace('.vcf', '')] + " ./" + file.getBaseName() + ".vcf" }.filter { it != null }.collect()
        knownSitesArgs_ch.view()

        filtered_vcf_ch = variantRecalibrator(final_vcf_ch, knownSitesArgs_ch, indexed_genome_ch.collect(), qsrc_vcf_ch.collect())
    } else {
        filtered_vcf_ch = filterVCF(final_vcf_ch, indexed_genome_ch.collect())
    }

    // Conditionally run identityAnalysis if identity_analysis is true
    if (params.identity_analysis) {

        // Create psam_info_ch and collect the sample ID and sex info into a single channel
        psam_info_ch = Channel
            .fromPath(params.samplesheet)
            .splitCsv(sep: '\t')
            .map { row -> tuple(row[0], row[3]) }

        // Initialize a variable to hold the combined PSAM content, starting with the header
        def combined_psam_content = new StringBuilder("#IID\tSID\tPAT\tMAT\tSEX\n")

        // Create a channel that processes sample information and appends it to the combined PSAM content
        psam_file_ch = psam_info_ch.map { sample_info ->
            def sample_id = sample_info[0]
            def sex = sample_info[1]

            // Convert empty sex values to 'NA' for unknown
            if (!sex) { sex = "NA" }

            // Generate PSAM content for this sample and strip any newlines before appending
            def sample_line = "${sample_id}\t${sample_id}\t0\t0\t${sex}".stripIndent().trim()
            combined_psam_content.append(sample_line + "\n")
        }

        // Save the combined content to a single .psam file and return the file path through the channel
        psam_file_ch.subscribe {
            def combined_psam_file = new File("/tmp/combined_samples.psam")
            combined_psam_file.text = combined_psam_content.toString()

            // Pass the file itself to the channel
            return combined_psam_file
        }

        // Now pass the psam_info_ch to the identityAnalysis process
        identity_analysis_ch = identityAnalysis(filtered_vcf_ch, psam_file_ch)
        identity_analysis_ch.view()
    }
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
