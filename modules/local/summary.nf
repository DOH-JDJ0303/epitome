process SUMMARY {
    label 'process_low'
    tag "${prefix}"

    input:
    tuple val(taxon), val(segment), path(qc_json), path(clusters_json), path(meta_csv), val(method)

    output:
    tuple val(taxon), val(segment), path("${prefix}.jsonl.gz"),   emit: json
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}.${method}"
    script = "epitome_summary.py"
    """
    # run script
    ${script} \\
        --taxon "${taxon}" \\
        --segment "${segment}" \\
        --method ${method} \\
        --qc ${qc_json} \\
        --clusters ${clusters_json} \\
        --meta ${meta_csv}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${script}: "\$(${script} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
