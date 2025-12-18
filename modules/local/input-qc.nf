process INPUT_QC {
    tag "${prefix}"
    label 'process_high'
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(sequences), path(metadata), path(exclusions)

    output:
    tuple val(taxon), val(segment), path("${prefix}.qc.fa.gz"),    emit: seqs, optional: true
    tuple val(taxon), val(segment), path("${prefix}.qc.jsonl.gz"), emit: summary
    path "versions.yml",                                           emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    tool = "epitome_qc.py"
    """
    ${tool} \\
        --taxon "${taxon}" \\
        --segment "${segment}" \\
        --z_threshold ${params.z_threshold} \\
        --amb_threshold ${params.amb_threshold} \\
        --fasta "${sequences}" \\
        ${metadata ? "--metadata ${metadata}" : ''}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
