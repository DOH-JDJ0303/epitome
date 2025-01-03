process EXPORT_QC_SEQS {
    tag "${prefix}"
    label 'process_low'
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), val(sequences)

    output:
    tuple val(taxon), val(segment), path("${prefix}.all.fa"), path("${prefix}.top.fa"), path("${prefix}.remainder.fa"), env(remainder_count), emit: seqs
    

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon}-${segment}"
    """
    # output sequences
    ## shuffle and number sequences
    cat ${sequences} | shuf | awk -v OFS='\t' '{print ">"\$1, \$2}' > shufd
    ## get seq count
    n_clean=\$(cat shufd | wc -l)
    ## parition sequences - this is necessary for large datasets
    cat shufd | tr '\t' '\n' > ${prefix}.all.fa
    cat shufd | sed -n "1,${params.max_cluster}p" | tr '\t' '\n' > ${prefix}.top.fa
    cat shufd | sed -n "${params.max_cluster+1},\\\$p" | tr '\t' '\n' > ${prefix}.remainder.fa
    remainder_count=\$(cat ${prefix}.remainder.fa | wc -l)
    """
}
