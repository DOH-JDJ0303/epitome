process SEQTK_SUBSEQ {
    tag "${prefix}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'quay.io/biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(taxa), val(segment), val(cluster), val(contigs), path(sequences), val(seq_count)

    output:
    tuple val(taxa), val(segment), val(cluster), path("${prefix}.fa"), emit: sequences
    path "versions.yml",                                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix = "${taxa}-${segment}-${cluster}"
    """
    echo "${contigs}" | tr -d '[],' | tr ' ' '\n' > contigs.txt
    seqtk \\
        subseq \\
        $args \\
        $sequences \\
        contigs.txt > ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = "fa"
    if ("$sequences" ==~ /.+\.fq|.+\.fq.gz|.+\.fastq|.+\.fastq.gz/) {
        ext = "fq"
    }
    """
    echo "" | gzip > ${sequences}${prefix}.${ext}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
