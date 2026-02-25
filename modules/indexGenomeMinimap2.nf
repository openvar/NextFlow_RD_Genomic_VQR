/*
 * Index the reference genome using Minimap2
 */

process indexGenomeMinimap2 {

    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }

    container 'quay.io/biocontainers/minimap2:2.28--h5bf99c6_0'

    tag "${referenceGenome.simpleName}"

    input:
    path referenceGenome

    output:
    path "${referenceGenome.simpleName}.mmi"

    script:
    """
    echo "Indexing reference genome ${referenceGenome.simpleName} with Minimap2"
    minimap2 -d ${referenceGenome.simpleName}.mmi ${referenceGenome}
    echo "Minimap2 genome indexing complete"
    """
}