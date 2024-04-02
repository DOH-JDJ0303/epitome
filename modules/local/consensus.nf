process CONSENSUS {
    tag "${prefix}"
    label 'process_high'
    container 'docker.io/jdj0303/bigbacter-base:1.0.0'

    input:
    tuple val(taxa), val(segment), val(cluster), path(aln)

    output:
    tuple val(taxa), val(segment), val(cluster), path("${prefix}.fa"), env(length), emit: fa
    path "${prefix}_length.csv",                                                    emit: len
    path "versions.yml",                                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}-${cluster}"

    """
    # run script
    consensus.sh "${prefix}" ${aln} "${task.cpus}"
    # collect consensus size info
    length=\$(cat ${prefix}.fa | grep -v '>' | tr -d '\n\t ' | wc -c)
    echo "${prefix},\${length}" > ${prefix}_length.csv

    # something about the normal way of getting version info messes with the creations of .command.env
    echo -e "\"${task.process}\":\n    consensus: \$(consensus.sh version)" > versions.yml
    """
}
