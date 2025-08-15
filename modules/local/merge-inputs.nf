process MERGE_INPUTS {
    tag "${prefix}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(assembly), path(metadata)

    output:
    tuple val(taxon), val(segment), path("${prefix}.assembly.fa.gz"), path("${prefix}.metadata.jsonl.gz"), emit: merged
    path "versions.yml",                                                                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    script = "epitome_metadata.py"
    """
    # combine assembly files
    gzip ${assembly.join(' ')} || true
    cat *.gz > ${prefix}.assembly.fa.gz

    # combine metadata
    ${script} \\
        --taxon "${taxon}" \\
        --segment "${segment}" \\
         ${metadata}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${script}: "\$(${script} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
