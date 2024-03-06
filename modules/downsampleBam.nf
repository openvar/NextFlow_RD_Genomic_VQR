/*
 * Downsample the BAM files
 */
process downsampleBam {

    label 'process_medium'
    container 'variantvalidator/indexgenome:1.1.0'

    tag "bamfile: $bamFile; fraction: $fraction"

    input:
    tuple val(fraction), file(bamFile)

    output:
    tuple val(bamFile.baseName), file("*_downsampled_*.bam")

    script:
    """
    echo "Running Downsampling"

    picard DownsampleSam \\
    I=${bamFile} \\
    O=${bamFile.baseName}_downsampled_${fraction}.bam \\
    P=${fraction} \\
    QUIET=true \\
    VALIDATION_STRINGENCY=SILENT > ${bamFile.baseName}_downsampled_${fraction}.bam.log
  
    echo "Downsampling complete"
    """
}
