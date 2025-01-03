process SUMMARY {
    label 'process_low'
    tag "${prefix}"

    input:
    tuple val(taxa), val(segment), path(ref_clusters), path(all_clusters), path(clean), path(raw), path(metadata), path(ani)

    output:
    tuple val(taxa), val(segment), path("*.csv"), emit: summary
    tuple val(taxa), val(segment), path("*.jpg"), emit: plots, optional: true
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    # convert reference and inputs into tab-separated format
    zcat ${raw} | cut -f 1 -d ' ' | sed 's/>.*\$/@&@/g' | tr -d '\n' | tr '@' '\n' | tail -n +2 | tr -d '>' | paste - - | awk -v OFS='\t' '{print \$1,toupper(\$2)}' > raw.txt
    cat ${clean} | sed 's/>.*\$/@&@/g' | tr -d '\n' | tr '@' '\n' | tail -n +2 | tr -d '>' | paste - - | awk -v OFS='\t' '{print \$1,toupper(\$2)}' > clean.txt

    # run script
    summary.R "${ref_clusters}" "${all_clusters}" raw.txt clean.txt "${metadata}" "${ani}" "${prefix}"

    # something about the normal way of getting version info messes with the creations of .command.env
    echo -e "\\"${task.process}\\":\\n    summary.R: \$(summary.R version)" > versions.yml
    """
}
