process INPUT_QC {
    tag "${prefix}"
    label 'process_low'
    container 'docker.io/jdj0303/bigbacter-base:1.0.0'
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(sequences)

    output:
    tuple val(taxa), val(segment), path("${prefix}.clean.fa"), env(seq_count), emit: assemblies
    path "${prefix}-qc-summary.csv",             emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"

    """
    seqs=${sequences}
    # check if file is compressed
    if [[ ${sequences} =~ \.gz\$ ]];
    then
        gzip -d ${sequences}
        seqs=\${seqs%.gz}
    fi
    # filter sequences
    input-qc.sh \${seqs}  ${prefix}
    # set sequence count
    seq_count=\$(cat ${prefix}-qc-summary.csv | cut -f 4 -d ',' | grep -v 'filter3' | tr -d '\t\r\n ')
    """
}
