process SUMMARY {
    label 'process_low'
    tag "${prefix}"

    input:
    tuple val(taxon), val(segment), path(qc_json), path(clusters_json), path(condensed_json), path(meta_csv)

    output:
    tuple val(taxon), val(segment), path("${prefix}.summary.json"),   emit: json
    tuple val(taxon), val(segment), path("${prefix}.summary.csv"), emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    """
    # run script
    epitome-summary.py \\
        --qc ${qc_json} \\
        --clusters ${clusters_json} \\
        --condensed ${condensed_json} \\
        --meta ${meta_csv}

    mv summary.json ${prefix}.summary.json
    mv summary.csv ${prefix}.summary.csv

    # something about the normal way of getting version info messes with the creations of .command.env
    echo -e "\\"${task.process}\\":\\n    epitome-summary.py: \$(epitome-summary.py --version)" > versions.yml
    """
}
