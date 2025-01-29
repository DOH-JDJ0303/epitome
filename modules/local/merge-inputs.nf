process MERGE_INPUTS {
    tag "${prefix}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(assembly), path(metadata)

    output:
    tuple val(taxon), val(segment), path("${prefix}.assembly.fa.gz"), path("${prefix}.metadata.csv"), emit: merged

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    """
    # combine assembly files
    gzip ${assembly.join(' ')} || true
    cat *.gz > ${prefix}.assembly.fa.gz
    # combine metadata
    merge-tables.R && mv merged.csv ${prefix}.metadata.csv
    """
}
