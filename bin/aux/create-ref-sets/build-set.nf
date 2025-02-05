#!/usr/bin/env nextflow

def currentDate = new Date().format('yyyy-MM-dd')

params.name  = params.name ? params.name : "EPITOME_${currentDate}"

workflow {
    // Gather taxon-segment options
    Channel
        .fromPath(params.input)
        .splitCsv(header: true, quote: '"')
        .map{ getChildren(it.run, ['pipeline_info']) }
        .flatten()
        .map{ getChildren(it), ['ncbi-data'] }
        .flatten()
        .set{ ch_input }
    
    // Pull reference genomes
    PULL_REFS (
        ch_input.map{ it.resolve('consensus') }.unique()
    )

    // Gather taxon-segment summaries
    ch_input
      .map{ it.list().toList().findAll{ f -> f =~ '.summary_full.csv' }.collect{ f -> it.resolve( f ) } }
      .flatten()
      .set{ ch_summary }

    // Merge summaries
    MERGE_SUMMARY (
        ch_summary.collect()
    )

    // Manage excluded sequences
    if(params.exclusions){
        MANAGE_EXCLUSIONS (
            file(params.exclusions)
        )
    }

    // Compress to tar.gz
    COMPRESS (
        MERGE_SUMMARY
            .out
            .result
            .map{ csv, md, acc -> [ csv, acc ] },
        PULL_REFS.out.result.collect(),
        params.exclusions ? MANAGE_EXCLUSIONS.out.result : []
    )

}

// Functions
def getChildren(dir, exclusions){
    def children = file(dir, checkIfExists: true).list().toList()
    children = children - exclusions
    return children.collect{ it -> file(dir).resolve( it ) }
}

// Processes
process MERGE_SUMMARY {
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.md"

    input:
    path summaries

    output:
    tuple path("refsheet.csv"), path("*.md"), path("accessions/*"), emit: result

    script:
    """
    merge-summary.R
    """
}
process PULL_REFS {
    tag "${consenus_dir}"

    input:
    path consenus_dir

    output:
    path "references/*", includeInputs: true, emit: result

    script:
    """
    mv consensus references
    """
}
process MANAGE_EXCLUSIONS {
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.md"

    input:
    path exclusions

    output:
    path "*", includeInputs: true, emit: result

    script:
    """
    manage-exclusions.R
    """
}
process COMPRESS {
    publishDir "${params.outdir}", mode: 'copy'
    stageInMode 'copy'

    input:
    tuple path(csv), path(acc, stageAs: 'accessions/*')
    path ref, stageAs: 'references/*'
    path exclusions

    output:
    path "*.tar.gz"

    script:
    """
    # compress the accessions separately
    tar cvzf accessions.tar.gz accessions/ && rm -r accessions
    # move everything to single directory & compress
    mkdir ${params.name} || true
    mv * ${params.name}/ || true
    tar cvzf ${params.name}.tar.gz ${params.name}/
    """
}