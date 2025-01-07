process MERGE_INPUTS {
    tag "${prefix}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(assembly), path(metadata)

    output:
    tuple val(taxa), val(segment), path("${prefix}.assembly.fa.gz"), path("${prefix}.metadata.csv"), emit: merged

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    # combine assembly files
    gzip ${assembly.join(' ')} || true
    cat *.gz > ${prefix}.assembly.fa.gz
    # combine metadata
    merge-tables.R && mv merged.csv ${prefix}.metadata.csv
    """
}
