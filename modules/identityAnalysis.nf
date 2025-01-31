process identityAnalysis {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/identity_analysis:0.0.1'

    tag "$vcfFile"

    publishDir("$params.outdir/IDENTITY_ANALYSIS", mode: "copy")

    input:
    tuple val(sample_id), file(vcfFile)              // Sample ID and VCF file
    file(psam_file)                                  // The pre-generated PSAM file

    output:
    path("relatedness_*")

    script:
    """
    echo "Running Identity Analysis for Sample: ${vcfFile}"

    # Step 1: Compress the VCF file using bgzip
    bgzip -c ${vcfFile} > ${vcfFile}.gz

    # Step 2: Index the compressed VCF file using tabix
    tabix -p vcf ${vcfFile}.gz

    # Step 3: Filter out multiallelic variants from the compressed VCF
    bcftools view -m2 -M2 ${vcfFile}.gz -o filtered_${sample_id}.vcf

    # Step 4: Convert the filtered VCF to PLINK format (for genetic analysis)
    plink2 --vcf filtered_${sample_id}.vcf --make-bed --split-par b38 --allow-extra-chr --update-sex ${psam_file} --out samples_${sample_id}

    # Step 5: Perform relatedness analysis with Plink using KING method
    plink2 --bfile samples_${sample_id} --king-cutoff 0.088 --make-king-table --allow-extra-chr --out relatedness_${sample_id}

    echo "Identity Analysis for Sample: ${vcfFile} Complete"
    """
}
