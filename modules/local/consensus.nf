process CONSENSUS {
    tag "${prefix}"
    label 'process_high'
    container 'docker.io/jdj0303/bigbacter-base:1.0.0'

    input:
    tuple val(taxa), val(segment), val(cluster), path(aln)

    output:
    tuple val(taxa), val(segment), val(cluster), path("${prefix}.fa"), emit: fa
    path "${prefix}_length.csv",                                       emit: len

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}-${cluster}"

    """
    # run script
    consensus.sh "${prefix}" ${aln} "${task.cpus}"
    # collect consensus size info
    echo "${prefix},\$(cat ${prefix}.fa | grep -v '>' | tr -d '\n\t ' | wc -c)" > ${prefix}_length.csv
    """
}
