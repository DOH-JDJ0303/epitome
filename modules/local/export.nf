process EXPORT {
    label 'process_low'

    input:
    path sheet
    val timestamp

    output:
    path "*.csv",        emit: samplesheet 
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # get outdir path - this is all possible in Groovy but I found it cleaner here
    path="${params.outdir}"
    path="\${path%/}"
    echo 'taxa,segment,assembly' > ${timestamp}-epitome.csv
    cat $sheet | awk -v OFS=',' -v p=\$path '{print \$1,\$2,p"/"\$1"/"\$2"/consensus/"\$3}' >> ${timestamp}-epitome.csv

    # version
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        export: 1.0
    END_VERSIONS
    """
}
