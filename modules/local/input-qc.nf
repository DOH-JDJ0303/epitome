process INPUT_QC {
    tag "${prefix}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(sequences), path(exclusions)

    output:
    tuple val(taxon), val(segment), path("${prefix}.qc.fa.gz"), emit: seqs
    tuple val(taxon), val(segment), path("${prefix}.qc.json"),  emit: summary
    path "versions.yml",                                        emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    """
    # gather metrics
    epitome-qc.py \\
        --taxon "${taxon}" \\
        --segment "${segment}" \\
        --amb_threshold ${params.amb_threshold} \\
        --len_threshold ${params.len_threshold} \\
        --fasta "${sequences}"

    gzip *.fa

    # odd stuff going on with versioning
    echo -e "\\"${task.process}\\":\\n    epitome-qc.py: \$(epitome-qc.py --version)" > versions.yml
    """
}
