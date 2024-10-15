process identityAnalysis {

    label 'process_medium'
    container 'your_custom_docker_with_plink_and_mtoolbox'

    tag "$vcfFile"

    publishDir("$params.outdir/identity_analysis", mode: "copy")

    input:
    tuple val(sample_id), file(vcfFile)
    file(mitoRef) optional true  // Optional mitochondrial reference file

    output:
    tuple val(sample_id), file("*.bed"), file("*.bim"), file("*.fam"),
          file("relatedness.genome"), file("roh_results.hom"),
          file("output_folder/*")

    script:
    """
    echo "Running Identity Analysis for Sample: ${vcfFile}"

    plink2 --vcf ${vcfFile} --make-bed --out samples
    plink2 --bfile samples --genome --out relatedness
    plink2 --bfile samples --homozyg --out roh_results

    # Extract mitochondrial variants
    bcftools view -r MT ${vcfFile} -o mtDNA_variants_sample1.vcf
    bcftools view -r MT sample2.vcf -o mtDNA_variants_sample2.vcf

    # Run MToolBox with rCRS reference
    MToolBox.sh -i mtDNA_variants_sample1.vcf mtDNA_variants_sample2.vcf -o output_folder --ref rCRS

    echo "Identity Analysis for Sample: ${vcfFile} Complete"
    """
}
