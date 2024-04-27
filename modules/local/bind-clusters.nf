process BIND_CLUSTERS {
    tag "${prefix}"
    label "process_low"
    container 'docker.io/johnjare/spree:1.0'

    input:
    tuple val(taxa), val(segment), path(clusters)

    output:
    path "*.csv",        emit: results


    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${taxa}-${segment}"
    """
    header="seq,taxa,segment,cluster"
    if [ -f "${prefix}-looseends.csv " ]
    then
        max_n=\$(cat *main.csv | cut -f 4 -d ',' | sort -rn)
        mv ${prefix}-looseends.csv tmp && cat tmp | grep -v "\$header" | tr ',' '\t' | awk -v max=\${max_n} -v OFS=',' '{print $1,$2,$3,$4+max}' > ${prefix}-looseends.csv
        rm tmp
    fi
    echo \$header > ${prefix}-clusters.csv
    cat ${clusters} | grep -v "\$header" >> ${prefix}-clusters.csv
    """
}
