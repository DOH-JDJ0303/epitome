process CLUSTER {
    tag "${prefix}"
    label 'process_low'
    container 'docker.io/johnjare/spree:1.0'

    input:
    tuple val(taxa), val(segment), path(dist)

    output:
    path "*.csv", emit: results
    path "*.jpg", emit: plot

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    # run script
    cluster.R ${dist} "${taxa}" "${segment}" ${params.threshold}
    """
}
