process SUMMARY {
    label 'process_low'
    tag "${prefix}"

    input:
    tuple val(taxa), val(segment), path(input_qc), path(clusters), path(refs), path(metadata)

    output:
    tuple val(taxa), val(segment), path("${prefix}.summary.csv"), emit: summary
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    # run script
    summary.R "${input_qc}" "${clusters}" "${refs}" "${metadata}"
    mv summary.csv ${prefix}.summary.csv

    # something about the normal way of getting version info messes with the creations of .command.env
    echo -e "\\"${task.process}\\":\\n    summary.R: \$(summary.R version)" > versions.yml
    """
}
