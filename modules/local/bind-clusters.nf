process BIND_CLUSTERS {
    tag "${prefix}"
    label "process_low"
    container 'docker.io/johnjare/spree:1.0'

    input:
    tuple val(taxa), val(segment), path(clusters)

    output:
    path "*.csv",        emit: results


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    header="seq,taxa,segment,cluster"
    echo \$header > ${prefix}-clusters.csv
    cat ${clusters} | grep -v "\$header" >> ${prefix}-clusters.csv
    """
}
