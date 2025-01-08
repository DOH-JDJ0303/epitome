process MERGE_CLUSTERS {
    tag "${prefix}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(main), path(assigned), path(looseends)

    output:
    tuple val(taxa), val(segment), path("${prefix}.clusters.csv"), emit: merged

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    merge-clusters.R "${main}" "${assigned}" "${looseends}"
    mv clusters.csv ${prefix}.clusters.csv
    """
}
