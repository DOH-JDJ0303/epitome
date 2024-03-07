process INPUT_QC {
    label 'process_low'
    container 'docker.io/jdj0303/bigbacter-base:1.0.0'

    input:
    path assemblies

    output:
    path "*.fa",                 emit: assemblies
    path "input-qc-summary.csv", emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # filter sequences
    input-qc.sh ${assemblies} #&& rm ${assemblies}
    """
}
