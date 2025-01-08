process CLUSTER {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(dist)
    val iteration

    output:
    tuple val(taxon), val(segment), path("*.csv"), emit: results
    path "*.jpg",        emit: plot
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    """
    zcat ${dist} | cut -f 1,2,3 > dists.txt
    # run script
    cluster.R dists.txt "${taxon}" "${segment}" "${iteration}" ${params.dist_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cluster.R: \$(cluster.R version)
    END_VERSIONS
    """
}