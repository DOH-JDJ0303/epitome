process CLUSTER {
    label 'process_low'
    container 'docker.io/johnjare/spree:1.0'

    input:
    path dist

    output:
    path "clusters.csv", emit: results
    path "clusters.jpg", emit: plot

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # run script
    cluster.R ${dist} ${params.threshold}
    """
}
