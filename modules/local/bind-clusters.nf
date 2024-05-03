process BIND_CLUSTERS {
    tag "${prefix}"
    label "process_low"
    container 'docker.io/johnjare/spree:1.0'

    input:
    tuple val(taxa), val(segment), path(clusters)

    output:
    tuple val(taxa), val(segment), path("*.csv"), emit: results


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    # define header
    header="seq,taxa,segment,cluster"
    # make loose-end cluster numbers continue from the firt (main) clustering.
    if [ -f "${prefix}-looseends.csv" ]
    then
        max_n=\$(cat *main.csv | cut -f 4 -d ',' | sort -rn)
        mv ${prefix}-looseends.csv tmp && cat tmp | grep -v "\$header" | tr ',' '\t' | awk -v m="\${max_n}" -v OFS=',' '{print \$1,\$2,\$3,\$4+m}' > ${prefix}-looseends.csv
        rm tmp
    fi
    # create combined file
    echo \$header > ${prefix}-clusters.csv
    cat ${clusters} | grep -v "\$header" >> ${prefix}-clusters.csv
    """
}
