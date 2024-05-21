process SUMMARY {
    label 'process_low'

    input:
    path summary
    path ani_ava
    path ani_seeds
    path seeds
    val timestamp


    output:
    path "*.csv",        emit: summary
    path "*.jpg",        emit: plots, optional: true
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # run script
    summary.R "${summary}" "${ani_ava}" "${ani_seeds}" "${seeds}" "${timestamp}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        summary: \$(summary.R version)
    END_VERSIONS
    """
}
