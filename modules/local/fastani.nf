process FASTANI_AVA {
    tag "${prefix}"
    label 'process_medium'

    conda "bioconda::fastani=1.32"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.32--he1c1bb9_0' :
        'biocontainers/fastani:1.32--he1c1bb9_0' }"

    input:
    tuple val(taxa), val(segment), path(assemblies), val(fraglen)

    output:
    tuple val(taxa), val(segment), path("${prefix}_ani.txt"), emit: ani
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${taxa}-${segment}"
    """
    # create list of assemblies
    echo "${assemblies.join('\n')}" > seqs.txt

    # run FastANI
    fastANI \\
        $args \\
        -t ${task.cpus} \\
        --fragLen $fraglen \\
        --ql seqs.txt \\
        --rl seqs.txt \\
        -o ${prefix}_ani.txt

    # odd stuff going on with versioning
    echo -e "\\"${task.process}\\":\\n    fastani: \$(fastANI --version 2>&1 | sed 's/version//;')" > versions.yml
    """
}

process FASTANI_SEEDS {
    label 'process_medium'

    conda "bioconda::fastani=1.32"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.32--he1c1bb9_0' :
        'biocontainers/fastani:1.32--he1c1bb9_0' }"

    input:
    path consensus, stageAs: 'consensus/*'
    path seeds, stageAs: 'seeds/*'

    output:
    path "seeds_ani.txt", emit: ani
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # create list of assemblies
    ls -d consensus/* > qfile.txt
    ls -d seeds/* > rfile.txt

    # get minimum contig length
    fragLen=\$(cat consensus/* seeds/* | sed 's/>.*\$/@&@/g' | tr -d '\n' | tr '@' '\n' | grep -v '>' | tail -n +2 | awk '{print length(\$0)-1}' | sort -n | sed -n 1p)

    # run FastANI
    fastANI \\
        $args \\
        -t ${task.cpus} \\
        --fragLen \${fragLen} \\
        --ql qfile.txt \\
        --rl rfile.txt \\
        -o seeds_ani.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
    """
}
