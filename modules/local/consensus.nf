process CONSENSUS {
    label 'process_high'
    container 'docker.io/jdj0303/bigbacter-base:1.0.0'

    input:
    tuple val(cluster), path(aln)

    output:
    path "${cluster}.fa", emit: fa

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # run script
    consensus.sh "${cluster}" ${aln}
    """
}
