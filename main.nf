// Use newest nextflow dsl
nextflow.enable.dsl = 2

// Print pipeline configuration
log.info """\
    ============================================
          DNASeq Pipeline Configuration
    ============================================
    samplesheet     : ${params.samplesheet}
    genome          : ${params.genome_file}
    genome index    : ${params.genome_index_files}
    index genome    : ${params.index_genome}
    qsr truth vcfs  : ${params.qsrVcfs}
    output directory: ${params.outdir}
    fastqc          : ${params.fastqc}
    aligner         : ${params.aligner}
    variant caller  : ${params.variant_caller}
    bqsr            : ${params.bqsr}
    degraded_dna    : ${params.degraded_dna}
    variant_recalibration: ${params.variant_recalibration}
    identity_analysis: ${params.identity_analysis}
    ============================================
""".stripIndent()

// Conditionally include modules
if (params.index_genome) {
    include { indexGenome } from './modules/indexGenome'
}
if (params.fastqc) {
    include { FASTQC } from './modules/FASTQC'
}
include { sortBam } from './modules/sortBam'
include { markDuplicates } from './modules/markDuplicates'
include { indexBam } from './modules/indexBam'
if (params.bqsr) {
    include { baseRecalibrator } from './modules/BQSR'
}
include { combineGVCFs } from './modules/processGVCFs'
include { genotypeGVCFs } from './modules/processGVCFs'
if (params.variant_recalibration) {
    include { variantRecalibrator } from './modules/variantRecalibrator'
} else {
    include { filterVCF } from './modules/filterVCF'
}
if (params.identity_analysis) {
    include { identityAnalysis } from './modules/identityAnalysis'
}
if (params.aligner == 'bwa-mem') {
    include { alignReadsBwaMem } from './modules/alignReadsBwaMem'
} else if (params.aligner == 'bwa-aln') {
    include { alignReadsBwaAln } from './modules/alignReadsBwaAln'
} else {
    error "Unsupported aligner: ${params.aligner}. Please specify 'bwa-mem' or 'bwa-aln'."
}
if (params.variant_caller == 'haplotype-caller') {
    include { haplotypeCaller } from './modules/haplotypeCaller'
} else {
    error "Unsupported variant caller: ${params.variant_caller}. Please specify 'haplotype-caller'."
}

if (params.degraded_dna) {
    include { mapDamage2 } from './modules/mapDamage'
    include { indexMapDamageBam } from './modules/indexBam'
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
        .map { row ->
            if (row.size() == 4) {
                tuple(row[0], [row[1], row[2]])
            } else if (row.size() == 3) {
                tuple(row[0], [row[1]])
            } else {
                error "Unexpected row format in samplesheet: $row"
            }
        }
    read_pairs_ch.view()

    // Run FASTQC on read pairs
    if (params.fastqc) {
        FASTQC(read_pairs_ch)
    }

    // Align reads to the indexed genome
    if (params.aligner == 'bwa-mem') {
        align_ch = alignReadsBwaMem(read_pairs_ch, indexed_genome_ch.collect())
    } else if (params.aligner == 'bwa-aln') {
        align_ch = alignReadsBwaAln(read_pairs_ch, indexed_genome_ch.collect())
    }

    // Sort BAM files
    sort_ch = sortBam(align_ch)

    // Mark duplicates in BAM files
    mark_ch = markDuplicates(sort_ch)

    // Index the BAM files and collect the output channel
    indexed_bam_ch = indexBam(mark_ch)

    // Conditionally run mapDamage if degraded_dna parameter is set
    if (params.degraded_dna) {
        // Run mapDamage2 process only if degraded_dna is true
        pre_mapDamage_ch = mapDamage2(indexed_bam_ch, indexed_genome_ch.collect())
        mapDamage_ch = indexMapDamageBam(pre_mapDamage_ch)
    } else {
        // If degraded_dna is not true, just pass through the sorted BAM files
        mapDamage_ch = indexed_bam_ch
    }

    // Create a channel from qsrVcfs
    knownSites_ch = Channel.fromPath(params.qsrVcfs)
        .filter { file -> file.getName().endsWith('.vcf.gz.tbi') || file.getName().endsWith('.vcf.idx') }
        .map { file -> "--known-sites " + file.getBaseName() }
        .collect()

    if (params.bqsr) {
        // Run BQSR on indexed BAM files
        bqsr_ch = baseRecalibrator(mapDamage_ch, knownSites_ch, indexed_genome_ch.collect(), qsrc_vcf_ch.collect())

    } else {
        // If BQSR is skipped, just pass through the mapDamage_ch channel
        bqsr_ch = mapDamage_ch
    }

    // Run HaplotypeCaller on BQSR files
    if (params.variant_caller == "haplotype-caller") {
        gvcf_ch = haplotypeCaller(bqsr_ch, indexed_genome_ch.collect()).collect()
    }

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
            'Homo_sapiens_assembly38.known_indels': 'known=true,training=false,truth=false,prior=15.0',  // High-priority known indels, not used for training
            'hapmap_3.3.hg38': 'known=false,training=false,truth=true,prior=15.0',  // Good for truth, not training
            '1000G_omni2.5.hg38': 'known=false,training=true,truth=false,prior=12.0',  // Omni SNPs, used for training
            '1000G_phase1.snps.high_confidence.hg38': 'known=true,training=true,truth=true,prior=10.0',  // High confidence SNPs, both for training and truth
            'Homo_sapiens_assembly38.dbsnp138': 'known=true,training=false,truth=false,prior=2.0',  // dbSNP, known but not for training
            'Mills_and_1000G_gold_standard.indels.hg38': 'known=true,training=true,truth=true,prior=12.0'  // Gold standard indels, good for truth (indels)
        ]

        // Generate --resource arguments
        knownSitesArgs_ch = Channel
            .fromPath(params.qsrVcfs)
            .filter { file -> file.getName().endsWith('.vcf.gz') || file.getName().endsWith('.vcf') }
            .map { file ->
                def baseName = file.getName().endsWith('.vcf.gz')
                    ? file.getName().replace('.vcf.gz', '')
                    : file.getName().replace('.vcf', '') // Correctly remove both .vcf.gz and .vcf
                def resourceArgs = resourceOptions.get(baseName) ?: "" // Get attributes from resourceOptions
                return "--resource:${baseName},${resourceArgs} ${file}" // Construct the full argument with proper formatting
            }
            .collect()
        filtered_vcf_ch = variantRecalibrator(final_vcf_ch, knownSitesArgs_ch, indexed_genome_ch.collect(), qsrc_vcf_ch.collect())
    } else {
        filtered_vcf_ch = filterVCF(final_vcf_ch, indexed_genome_ch.collect())
    }

    // Conditionally run identityAnalysis if identity_analysis is true
    if (params.identity_analysis) {

        //Create psam_info_ch and collect the sample ID and sex info into a single channel
        psam_info_ch = Channel
            .fromPath(params.samplesheet)
            .splitCsv(sep: '\t')
            .map { row ->
                if (row.size() == 4) {
                    tuple(row[0], row[3])  // Sample ID and sex info when 4 columns are present
                } else if (row.size() == 3) {
                    tuple(row[0], row[2])    // Sample ID and null for sex info when 3 columns are present
                } else {
                    error "Unexpected row format in samplesheet: $row"  // Handle unexpected formats
                }
            }


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
    }
}

workflow FASTQC_only {
    // Set channel to gather read_pairs
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row ->
            if (row.size() == 4) {
                tuple(row[0], [row[1], row[2]])
            } else if (row.size() == 3) {
                tuple(row[0], [row[1]])
            } else {
                error "Unexpected row format in samplesheet: $row"
            }
        }
    read_pairs_ch.view()

    if (params.fastqc) {
        FASTQC(read_pairs_ch)
    }
}

workflow.onComplete {
    log.info ( workflow.success ? "\nworkflow is done!\n" : "Oops .. something went wrong" )
}