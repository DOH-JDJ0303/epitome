process CONDENSE {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(dist), path(consensus), path(clusters)

    output:
    tuple val(taxa), val(segment), path("${prefix}.condensed.csv"), path("*.fa.gz", includeInputs: true), env(min_len), emit: results
    path "versions.yml",                                                                                                emit: versions


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
    for s in \$(cat ${prefix}.condensed.csv | tail -n +2 | cut -f 1 -d ',')
    do
        mv \${s}.fa.gz tmp/
    done
    rm *.fa.gz || true
    mv tmp/*.fa.gz ./ || true 
    rm -r tmp || true
    
    # get min seq length - for FastANI
    min_len=\$(cat ${prefix}.condensed.csv | cut -f 7 -d ',' | tail -n +2 | sort -n | sed -n 1p)

    # something about the normal way of getting version info messes with the creations of .command.env
    echo -e "\\"${task.process}\\":\\n    condense.R: \$(condense.R version)" > versions.yml
    """
}