process CONDENSE {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(dist), path(consensus), path(clusters)

    output:
    tuple val(taxa), val(segment), path("${prefix}.condensed.csv"), path("*.fa.gz", includeInputs: true), emit: results
    path "versions.yml",                                                                                  emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    # get seq lengths
    zcat ${consensus} | paste - - | tr -d '>' | sed 's/.fa//g' | awk -v OFS=',' '{print \$1,length(\$2)}' > lengths.csv

    # run script
    condense.R ${dist} lengths.csv ${clusters} "${taxa}" "${segment}" ${params.dist_threshold}

    # remove sequences that will not be retained
    mkdir tmp
    for s in \$(cat ${prefix}.condensed.csv | tr -d '"' | tail -n +2 | awk -v FS=',' '\$1 != "condensed" {print \$1}')
    do
        mv \${s}.fa.gz tmp/
    done
    rm *.fa.gz || true
    mv tmp/*.fa.gz ./ || true 
    rm -r tmp || true

    # something about the normal way of getting version info messes with the creations of .command.env
    echo -e "\\"${task.process}\\":\\n    condense.R: \$(condense.R version)" > versions.yml
    """
}