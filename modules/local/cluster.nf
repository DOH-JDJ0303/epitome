process CLUSTER {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(seqs)

    output:
    tuple val(taxon), val(segment), path("*.csv"), emit: results
    path "*.png",        emit: plot, optional: true
    // path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    """
    epitome-cluster.py \\
        --fasta ${seqs} \\
        --max_cluster ${params.max_cluster} \\
        --dist_threshold ${params.dist_threshold} \\
        --threads ${task.cpus}
    mv clusters.csv ${prefix}.clusters.csv
    """
}