process ASSIGN_REMAINDER {
    tag "${prefix}"
    label "process_high"
    container 'docker.io/johnjare/spree:1.0'
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(reps), path(mash), path(all)

    output:
    tuple val(taxa), val(segment), path("${prefix}-assigned.csv"), emit: assigned
    tuple val(taxa), val(segment), path("not-assigned.csv"), env(not_assigned_count), emit: not_assigned

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    assign-remainder.R ${reps} ${mash} ${all} ${taxa} ${segment}
    not_assigned_count=\$(cat not-assigned.csv | wc -l)
    """
}