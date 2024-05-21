process INPUT_QC {
    tag "${prefix}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(sequences), val(expected_length)

    output:
    tuple val(taxa), val(segment), path("${prefix}.top.fa"), path("${prefix}.remainder.fa"), env(remainder_count), emit: seqs
    tuple val(taxa), val(segment), path("${prefix}.all.fa"),                                                       emit: all
    path "${prefix}-qc-summary.csv",                                                                               emit: summary
    path "versions.yml",                                                                                           emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    seqs=${sequences}
    # check if file is compressed
    if [[ "${sequences}" == *".gz" ]];
    then
        gzip -d ${sequences}
        seqs=\${seqs%.gz}
    fi
    # filter sequences
    input-qc.sh \${seqs} ${prefix} "${expected_length}" "${params.len_threshold}" "${params.max_cluster}"

    cat ${prefix}.top.fa ${prefix}.remainder.fa > ${prefix}.all.fa

    # count remainder
    remainder_count=\$(cat ${prefix}.remainder.fa | paste - - | wc -l)

    # something about the normal way of getting version info messes with the creations of .command.env
    echo -e "\\"${task.process}\\":\\n    input-qc.sh: \$(input-qc.sh version)" > versions.yml
    """
}
