process SEQTK_SUBSEQ {
    tag "${prefix}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'quay.io/biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(taxa), val(segment), val(cluster), val(seqs), path(contigs)

    output:
    tuple val(taxa), val(segment), val(cluster), path("${prefix}.fa"), val("${seqs.size()}"), emit: sequences
    path "versions.yml",                                                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix = "${taxa}-${segment}-${cluster}"
    """
    seqtk \\
        subseq \\
        ${args} \\
        ${contigs} \\
        <(echo "${seqs.join('\n')}") > ${prefix}.fa

    # odd stuff going on with versioning
    echo -e "\\"${task.process}\\":\\n    seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')" > versions.yml
    """
}

process SEQTK_LOOSEENDS {
    tag "${prefix}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'quay.io/biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(taxa), val(segment), path(seq_list), path(sequences)

    output:
    tuple val(taxa), val(segment), path("${prefix}-looseends.fa"), emit: sequences
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix = "${taxa}-${segment}"
    """
    seqtk \\
        subseq \\
        $args \\
        $sequences \\
        $seq_list > ${prefix}-looseends.fa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
