process CLUSTER {
    tag "${prefix}"
    label "process_high"
    container 'docker.io/johnjare/spree:1.0'
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(dist), val(seq_count)

    output:
    path "*.csv",        emit: results
    path "*.jpg",        emit: plot
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    zcat ${dist} | cut -f 1,2,3 > dists.txt
    # run script
    cluster.R dists.txt "${taxa}" "${segment}" ${params.dist_threshold}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cluster: \$(cluster.R version)
        END_VERSIONS
    """
}

process CLUSTER_LARGE {
    tag "${prefix}"
    label "process_high_memory"
    container 'docker.io/johnjare/spree:1.0'
    errorStrategy 'retry'
    maxRetries 4

    input:
    tuple val(taxa), val(segment), path(dist), val(seq_count)

    output:
    path "*.csv",        emit: results
    path "*.jpg",        emit: plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    zcat ${dist} | cut -f 1,2,3 > dists.txt
    # run script
    cluster.R dists.txt "${taxa}" "${segment}" ${params.dist_threshold}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cluster: \$(cluster.sh version)
    END_VERSIONS
    """
}
