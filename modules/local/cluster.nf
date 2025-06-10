process CLUSTER {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(seqs)

    output:
    tuple val(taxon), val(segment), path("*.json"),  emit: json
    tuple val(taxon), val(segment), path("*.multi.fa.gz"), emit: multi, optional: true
    tuple val(taxon), val(segment), path("*.single.fa.gz"), emit: single, optional: true

    // path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    """
    epitome-cluster.py \\
        --fasta ${seqs} \\
        --taxon ${taxon} \\
        --segment ${segment} \\
        --max_cluster ${params.max_cluster} \\
        --dist ${params.dist_threshold} \\
        --window ${params.window_size} \\
        --ksize ${params.ksize} \\
        --scaled ${params.scaled}
    
    echo -e "\\"${task.process}\\":\\n    epitome-cluster.py: \$(epitome-cluster.py --version)" > versions.yml
    """
}