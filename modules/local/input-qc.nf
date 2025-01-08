process INPUT_QC {
    tag "${prefix}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(sequences)

    output:
    tuple val(taxon), val(segment), path("${prefix}.all.fa.gz"), path("${prefix}.top.fa.gz"), path("${prefix}.remainder.fa.gz"), env(STATUS), emit: seqs
    tuple val(taxon), val(segment), path("${prefix}.qc.csv"),                                                                                 emit: summary
    path "versions.yml",                                                                                                                      emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    """
    # decompress FASTA file, if needed
    gzip -d ${sequences} || true
    # gather metrics
    input-qc.R "${taxon}" "${segment}" "${sequences.name.replaceAll(/\.gz$/,'')}" "${params.amb_threshold}" "${params.len_threshold}" "${params.max_cluster}"
    gzip *.fa
    # check for top and remainder files
    if [ ! -e "${prefix}.top.fa.gz" ]
    then
        STATUS=null
        touch ${prefix}.top.fa.gz ${prefix}.remainder.fa.gz
    else
        STATUS=subsampled
    fi
    # odd stuff going on with versioning
    echo -e "\\"${task.process}\\":\\n    input-qc.R: \$(input-qc.R version)" > versions.yml
    """
}
