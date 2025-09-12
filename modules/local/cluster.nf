process CLUSTER {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(seqs)

    output:
    tuple val(taxon), val(segment), path("*.jsonl.gz"),     emit: json
    tuple val(taxon), val(segment), path("*.multi.fa.gz"),  emit: multi, optional: true
    tuple val(taxon), val(segment), path("*.single.fa.gz"), emit: single, optional: true
    path "versions.yml",                                    emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    tool = "epitome_cluster.py"
    """
    ${tool} \\
        --fasta ${seqs} \\
        --taxon "${taxon}" \\
        --segment "${segment}" \\
        --max_cluster ${params.max_cluster} \\
        --dist ${params.dist_threshold} \\
        --window_size ${params.window_size} \\
        --ksize ${params.ksize} \\
        --scaled ${params.scaled} \\
        ${params.centroid ? '--centroid' : ''}
    
    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}