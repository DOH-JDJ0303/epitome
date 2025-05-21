process CONSENSUS {
    tag "${prefix}"
    label 'process_high'
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), val(cluster), path(aln)

    output:
    tuple val(taxon), val(segment), path("${prefix}.fa.gz"), emit: fa
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}-${cluster}"
    """
    # run script
    consensus.sh "${cluster}" ${aln}
    # compress
    mv ${cluster}.fa ${prefix}.fa && gzip ${prefix}.fa
    # collect consensus size info
    length=\$(zcat ${prefix}.fa.gz | grep -v '>' | tr -d '\n\t ' | wc -c)
    echo "${prefix},\${length}" > ${prefix}_length.csv

    # something about the normal way of getting version info messes with the creations of .command.env
    echo -e "\\"${task.process}\\":\\n    consensus.sh: \$(consensus.sh version)" > versions.yml
    """
}
