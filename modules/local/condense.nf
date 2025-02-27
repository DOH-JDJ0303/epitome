process CONDENSE {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(seqs), path(clusters)

    output:
    tuple val(taxon), val(segment), path("${prefix}.condensed.csv"), path("*.fa.gz", includeInputs: true), emit: results
    tuple val(taxon), val(segment), path("condense,png"),                                                  emit: plot, optional: true
    // path "versions.yml",                                                                                emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    """
    # combine sequences into single file
    cat ${seqs} > seqs.fa.gz
    # run script
    epitome-condense.py \\
        --fasta seqs.fa.gz \\
        --clusters ${clusters} \\
        --dist_threshold ${params.dist_threshold}
    
    # rename output
    mv condensed.csv ${prefix}.condensed.csv
    # remove condensed sequences
    rm seqs.fa.gz
    mkdir condensed
    for s in \$(cat ${prefix}.condensed.csv | tr -d '"' | tail -n +2 | awk -v FS=',' '\$5 != "" {print \$4}' | uniq)
    do
        mv \${s}.fa.gz condensed/
    done
    """
}