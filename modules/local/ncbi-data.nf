process NCBI_DATA {
    tag "${taxon}"
    label 'process_low'

    input:
    tuple val(taxon)

    output:
    tuple val(taxon), path("ncbi_dataset/data/genomic.fna.gz"), emit: genomic
    tuple val(taxon), path("parsed/*.json"),                    emit: data_reports
    tuple val(taxon), path("ncbi-taxids.json"),                 emit: taxids
    tuple val(taxon), path("ncbi-subtype.csv"),                 emit: subtype
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args   = task.ext.args ?: ''
    """    
    # Download sequences & limited metadata
    datasets download virus genome taxon "${taxon}" ${args} && unzip ncbi_dataset.zip
    # compress sequence file
    gzip ncbi_dataset/data/genomic.fna
    # parse data into individual JSON files
    mkdir parsed
    awk '{print > "parsed/"NR".json"}' ncbi_dataset/data/data_report.jsonl
    
    # Download detailed taxonomy data
    datasets summary taxonomy taxon "${taxon}" --rank species > ncbi-taxids.json

    # Download additional details (those not included with NCBI Datasets)
    esearch -db nucleotide -query '${taxon}[Organism] AND "complete sequence"[Title]' | \
       efetch -format docsum -mode json | \
       jq -r '.result | .[] | [.accessionversion, .subtype, .subname]? | @csv' \
       > ncbi-subtype.csv

    # version info
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ncbi-datasets: \$(datasets --version | sed 's/.*: //g')
    END_VERSIONS
    """
}
