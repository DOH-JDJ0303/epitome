process MAFFT {
    tag "${taxon}-${segment}-${cluster}"
    label 'process_high'

    conda "bioconda::mafft=7.520"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6:425707898cf4f85051b77848be253b88f1d2298a-0':
        'quay.io/biocontainers/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6:425707898cf4f85051b77848be253b88f1d2298a-0' }"

    input:
    tuple val(taxon), val(segment), val(cluster), path(seqs)

    output:
    tuple val(taxon), val(segment), val(cluster), path("${prefix}.mafft.fa"), emit: fa
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    def args         = task.ext.args   ?: ''
    prefix = "${taxon.replaceAll(' ','_')}-${segment}-${cluster}"
    """
    # run mafft
    mafft \\
        ${args} \\
        --thread ${task.cpus} \\
        --auto \\
        ${seqs} \\
        > ${prefix}.mafft.fa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}