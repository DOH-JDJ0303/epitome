process MERGE_INPUTS {
    tag "${prefix}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(taxon), path(assembly), path(metadata), val(segmented)

    output:
    tuple val(taxon), path("*.fa.gz"),             emit: fa
    tuple val(taxon), path("*.metadata.jsonl.gz"), emit: meta
    tuple val(taxon), path("manifest.csv"),        emit: man
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    prefix = "${taxon.replaceAll(' ','_')}"
    tool = "epitome_inputs.py"
    """
    # combine metadata
    ${tool} \\
        --taxon "${taxon}" \\
        --assembly ${assembly} \\
        --metadata ${metadata} \\
        ${segmented ? '--segmented' : ''} \\
        ${args}


    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
