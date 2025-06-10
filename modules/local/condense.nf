process CONDENSE {
    tag "${prefix}"
    label "process_high"
    stageInMode 'copy'

    input:
    tuple val(taxon), val(segment), path(seqs), path(clusters)

    output:
    tuple val(taxon), val(segment), path("${prefix}.condensed.json"), path("*.fa.gz", includeInputs: true), emit: results
    tuple val(taxon), val(segment), path("condense.png"),                                                   emit: plot, optional: true
    path "versions.yml",                                                                                    emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxon.replaceAll(' ','_')}-${segment}"
    """
    mkdir condensed || mv condensed/* ./ || true

    # run script
    epitome-condense.py \\
        --taxon ${taxon} \\
        --segment ${segment} \\
        --clusters ${clusters} \\
        --dist ${params.dist_threshold} \\
        --window ${params.window_size} \\
        --ksize ${params.ksize} \\
        --scaled ${params.scaled} \\
        --fasta ${seqs}

    # move condensed (not published)
    for s in \$(cat *.condensed.json | tr -d '{":\t ' | grep -B 1 "condensed${prefix}" | paste - - - | cut -f 1 | uniq)
    do
        dis_file=\${s}.fa.gz
        echo "Discarding \$dis_file"
        mv \$dis_file condensed/
    done
    
    # update final sequence headers
    for file in *.fa.gz; do
    base=\$(basename "\$file" .fa.gz)
    zcat "\$file" | awk -v prefix="\$base" '{ if (\$0 ~ /^>/) { print ">" prefix } else { print } }' | gzip > tmp.fa.gz
    mv tmp.fa.gz \$file
    done

    # odd stuff going on with versioning
    echo -e "\\"${task.process}\\":\\n    epitome-condense.py: \$(epitome-condense.py --version)" > versions.yml
    """
}