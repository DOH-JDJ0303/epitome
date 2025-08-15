process CONSENSUS {
    tag "${prefix}"
    label 'process_high'
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), val(cluster), path(aln)

    output:
    tuple val(taxon), val(segment), path("${prefix}.fa.gz"), emit: fa
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}-${cluster}"
    script = "epitome_consensus.py"
    """
    # run script
    ${script} \\
        --prefix ${prefix} \\
        --aln ${aln}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${script}: "\$(${script} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
