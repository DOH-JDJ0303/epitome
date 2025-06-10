process NCBI_DATA {
    tag "${taxon}"
    label 'process_low'

    input:
    tuple val(taxon), val(segment_synonyms), val(seg_status)

    output:
    tuple val(taxon), path("genomic.fa.gz"),             emit: genomic
    tuple val(taxon), path("ncbi-meta.complete.csv"),     emit: complete
    tuple val(taxon), path("ncbi-meta.report-only.csv"),  emit: reportonly
    tuple val(taxon), path("ncbi-meta.subtype-only.csv"), emit: subtypeonly
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args   = task.ext.args ?: ''
    """    
    #---- NCBI Datasets Genome -----#
    # Download virus genomes 
    datasets download virus genome taxon "${taxon}" ${args} && unzip ncbi_dataset.zip
    # clean & compress sequence file
    cat ncbi_dataset/data/genomic.fna | awk '/^>/ {print \$1; next} {print}' | gzip > genomic.fa.gz
    
    #---- NCBI Datasets Taxon -----#
    # Download detailed taxonomy data
    datasets summary taxonomy taxon "${taxon}" --rank species > ncbi-taxids.json
    
    #---- NCBI Esearch -----#
    # Download additional details (those not included with NCBI Datasets)
    esearch -db nucleotide -query '${taxon}[Organism] AND "complete sequence"[Title]' | \
       efetch -format docsum -mode json \
       > ncbi-subtype.json
    
    #---- COMBINE ----#
    ncbi-data.py --segsyns "${segment_synonyms}" ${seg_status == 'TRUE' ? '--segmented' : ''}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        ncbi-datasets: \$(datasets --version | sed 's/.*: //g')
    END_VERSIONS
    """
}
