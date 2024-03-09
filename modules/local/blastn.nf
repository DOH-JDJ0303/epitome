process BLASTN {
    tag "${prefix}"
    label 'process_medium'

    conda "bioconda::blast=2.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.14.1--pl5321h6f7f691_0':
        'quay.io/biocontainers/blast:2.14.1--pl5321h6f7f691_0' }"

    input:
    tuple val(taxa), val(segment), val(cluster), path(consensus)

    output:
    path "${prefix}-blast.txt", emit: results
    path "versions.yml",        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = "${taxa}-${segment}"
    """
    # combine assemblies
    cat ${consensus} > all.fa
    # run Blastn
    blastn \\
        -num_threads ${task.cpus} \\
        -subject all.fa \\
        -query all.fa \\
        ${args} \\
        | awk -v OFS='\t' -v taxa="${taxa}" -v seg="${segment}" '{print taxa, seg, \$0}' \\
        > ${prefix}-blast.txt
    # clean up
    rm all.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
