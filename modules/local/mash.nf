process MASH_TOP {
    tag "${taxa}-${segment}"
    label 'process_high'
    conda "bioconda::mash=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'quay.io/biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(taxa), val(segment), path(sequences)

    output:
    tuple val(taxa), val(segment), path("${prefix}-dist.txt.gz"), emit: dist
    path "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${taxa}-${segment}"

    """
    # create mash sketch
    mash \\
        sketch \\
        $args \\
        -p $task.cpus \\
        -o sketch \\
        -i ${sequences}

    # calculate distance
    mash \\
        dist \\
        $args \\
        -p $task.cpus \\
        sketch.msh \\
        sketch.msh \\
        | gzip > ${prefix}-dist.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}

process MASH_REMAINDER {
    tag "${taxa}-${segment}"
    label 'process_high'
    conda "bioconda::mash=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'quay.io/biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(taxa), val(segment), path(top), path(remainder), path(clusters) 

    output:
    tuple val(taxa), val(segment), path("reps.csv"), path("remainder-mash.csv"), path("remainder-list.csv"), emit: results
    path "versions.yml",                                                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${taxa}-${segment}"

    """
    assign-remainder.sh ${top} ${remainder} ${clusters} ${params.dist_threshold} $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
