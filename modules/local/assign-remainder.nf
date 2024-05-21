process ASSIGN_REMAINDER {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(reps), path(mash), path(all)

    output:
    tuple val(taxa), val(segment), path("${prefix}-assigned.csv"),                    emit: assigned
    tuple val(taxa), val(segment), path("not-assigned.csv"), env(not_assigned_count), emit: not_assigned
    path "versions.yml",                                                              emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    # assign remainder
    assign-remainder.R ${reps} ${mash} ${all} ${taxa} ${segment}
    # get count of sequences not assigned
    not_assigned_count=\$(cat not-assigned.csv | wc -l)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assign-remainder.R: \$(assign-remainder.R version)
    END_VERSIONS
    """
}