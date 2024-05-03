process CONDENSE {
    tag "${prefix}"
    label "process_high"
    container 'docker.io/johnjare/spree:1.0'
    stageInMode 'copy'

    input:
    tuple val(taxa), val(segment), path(dist), path(consensus), path(clusters)

    output:
    tuple val(taxa), val(segment), path("${prefix}.condensed.csv"), path("*.fa", includeInputs: true), env(min_len), emit: results
    //path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    # get seq lengths
    cat ${consensus} | paste - - | tr -d '>' | sed 's/.fa//g' | awk -v OFS=',' '{print \$1,length(\$2)}' > lengths.csv

    # run script
    condense.R ${dist} lengths.csv ${clusters} "${taxa}" "${segment}" ${params.dist_threshold}

    # remove sequences that will not be retained
    mkdir tmp
    for s in \$(cat ${prefix}.condensed.csv | tail -n +2 | cut -f 1 -d ',')
    do
        mv \${s}.fa tmp/
    done
    rm *.fa && mv tmp/*.fa ./ && rm -r tmp
    
    # get min seq length - for FastANI
    cat ${prefix}.condensed.csv | cut -f 7 -d ',' | tail -n +2 | sort -n | sed -n 1p | awk '{print "min_len="\$1-1}' > .command.env

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cluster: \$(condense.R version)
        END_VERSIONS
    """
}