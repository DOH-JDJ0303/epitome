process CLUSTER {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(dist)
    val iteration

    output:
    tuple val(taxa), val(segment), path("*.csv"), emit: results
    path "*.jpg",        emit: plot
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    zcat ${dist} | cut -f 1,2,3 > dists.txt
    # run script
    cluster.R dists.txt "${taxa}" "${segment}" "${iteration}" ${params.dist_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cluster.R: \$(cluster.R version)
    END_VERSIONS
    """
}