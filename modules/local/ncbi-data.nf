process NCBI_DATA {
    tag "${taxon}"
    label 'process_low'

    input:
    tuple val(taxon), val(segment_synonyms), val(seg_status)

    output:
    tuple val(taxon), path("*.fa.gz"), emit: fa
    tuple val(taxon), path("*.json"),  emit: json
    tuple val(taxon), path("*.csv"),   emit: man
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args   = task.ext.args ?: ''
    script = "epitome_ncbi.py"
    """
    mkdir data/ || true

    #---- NCBI Datasets Genome -----#
    # Downloads genome assemblies and associated sample data.
    # Sample data is currently lacking subtype information (Downloaded below)

    datasets download virus genome taxon \\
        "${taxon}" \\
        ${args}
    
    unzip ncbi_dataset.zip
    mv ncbi_dataset/data/genomic.fna data/datasets-genome.fa
    mv ncbi_dataset/data/data_report.jsonl data/datasets-genome.jsonl
    
    #---- NCBI Datasets Taxon -----#
    # Used to identify the species-level classification

    datasets summary taxonomy taxon \\
        "${taxon}" \\
        --rank species \\
        > data/datasets-taxonomy.json
    
    #---- NCBI EDirect -----#
    # Gather subtype information (Not included in Datasets report)

    cat data/datasets-genome.fa \\
        | grep '>' \\
        | cut -f 1 -d ' ' \\
        | tr -d '>' \\
        > data/accessions.txt

    epost -db nuccore -input data/accessions.txt \\
        | efetch -format docsum -mode json \\
        > data/edirect.json
    
    #---- COMBINE ----#
    ${script} \\
        --taxon ${taxon.replace(' ', '_')} \\
        --segment_synonyms "${segment_synonyms}" \\
        ${seg_status == 'TRUE' ? '--segmented' : ''} \\
        --datasets_genome_fasta data/datasets-genome.fa \\
        --datasets_genome_json data/datasets-genome.jsonl \\
        --datasets_taxonomy data/datasets-taxonomy.json \\
        --edirect data/edirect.json

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets: \$(datasets --version | sed 's/.*: //g')
        ${script}: "\$(${script} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
