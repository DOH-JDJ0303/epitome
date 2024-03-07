process MASH {
    label 'process_medium'
    conda "bioconda::mash=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'biocontainers/mash:2.3--he348c14_1' }"

    input:
    path sequences

    output:
    path("dist.txt")   , emit: dist
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # create mash sketch
    mash \\
        sketch \\
        $args \\
        -p $task.cpus \\
        -o sketch \\
        $sequences


    # calculate distance
    mash \\
        dist \\
        $args \\
        -p $task.cpus \\
        sketch.msh \\
        sketch.msh \\
        > dist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.msh
    touch ${prefix}.mash_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
