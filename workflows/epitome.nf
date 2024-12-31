/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowRefmaker.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { TIMESTAMP                    } from '../modules/local/timestamp'
include { NCBI_DATA                    } from '../modules/local/ncbi-data'
include { INPUT_QC                     } from '../modules/local/input-qc'
include { MASH_TOP                     } from '../modules/local/mash'
include { CLUSTER as CLUSTER_MAIN      } from '../modules/local/cluster'
include { MASH_REMAINDER               } from '../modules/local/mash'
include { ASSIGN_REMAINDER             } from '../modules/local/assign-remainder'
include { SEQTK_LOOSEENDS              } from '../modules/local/seqtk_subseq'
include { MASH_TOP as MASH_LOOSEENDS   } from '../modules/local/mash'
include { CLUSTER as CLUSTER_LOOSEENDS } from '../modules/local/cluster'
include { SEQTK_SUBSEQ                 } from '../modules/local/seqtk_subseq'
include { MAFFT                        } from '../modules/local/mafft'
include { CONSENSUS                    } from '../modules/local/consensus'
include { MASH_TOP as MASH_CONDENSE    } from '../modules/local/mash'
include { CONDENSE                     } from '../modules/local/condense'
include { FASTANI_AVA                  } from '../modules/local/fastani'
include { SUMMARY                      } from '../modules/local/summary'    
include { EXPORT                       } from '../modules/local/export'    


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

def formatTaxonData (row) {
    def data = [:]
    row[2].each{ data[it.key] = it.value }

    def result = [ taxon: row[0], 
                   accession: data.accession, 
                   length: data.length, 
                   collectionDate: data.containsKey('isolate') ? ( data.isolate.containsKey('collectionDate') ? data.location.collectionDate : null ) : null, 
                   geographicRegion: data.containsKey('location') ? ( data.location.containsKey('geographicRegion') ? data.location.geographicRegion : null ) : null,
                   organismName_host: data.containsKey('host') ? ( data.host.containsKey('organismName') ? data.host.organismName : null ) : null,
                   organismName_virus: data.virus.organismName.split(' ').findAll{ ! it.contains('/') }.join(' '),
                   taxIds: data.virus.lineage.taxId
                   ]
    return result
}

def formatSubtypeData (row) {
    def target_keys = ['segment','subtype','genotype']
    def results     = [ accession: row[1][0] ]
    def keys        = row[1][1].split('\\|').toList()
    def values      = row[1][2].split('\\|').toList()
    if( keys.size() == values.size() ){
        values.eachWithIndex{ value, index -> results[ keys[index] ] = value }
    }
    target_keys.findAll{ ! results.keySet().contains(it) }.each{ results[it] = null }

    return  results.findAll { it -> ('accession'+target_keys).contains(it.key) }

}

def fixSegmentSynonyms(row){
    def taxon = row[0]
    def data  = row[1]
    def syns  = row[2]
    if(syns){
        data.segment = syns.findResults { key, value -> value.contains(data.segment) ? key : null }[0]
    }
    return [ taxon, data ]
}

workflow EPITOME {

    ch_versions = Channel.empty()

    /*
    =============================================================================================================================
        LOAD SAMPLESHEET & GET TIMESTAMP
    =============================================================================================================================
    */

    // MODULE: Get timestamp
    TIMESTAMP ()
        .timestamp
        .set{ch_timestamp}

    // Load samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(header:true)
        .map{ [ taxon:    it.taxon, 
                segment:  it.containsKey('segment') ? ( it.segment ? it.segment : 'wg') : 'wg',
                assembly: it.containsKey('assembly') ? ( it.assembly ? file(it.assembly, checkIfExists: true) : [] ) : [], 
                length:   it.containsKey('length') ? ( it.length ? it.length : null ) : null, 
                metadata: it.containsKey('metadata') ? ( it.metadata ? file(it.metadata, checkIfExists: true) : [] ) : [],
                segmentSynomyms: it.containsKey('segmentSynonyms') ? ( it.segmentSynonyms ? it.segmentSynonyms : null ) : null,
                ] 
                }
        .set{ ch_input }

    /*
    =============================================================================================================================
        GATHER TAXON DATA
    =============================================================================================================================
    */
    if(params.ncbi){
        NCBI_DATA(
            ch_input.map{ it.taxon }.unique()
        )
        
        NCBI_DATA
            .out
            .taxids
            .splitJson()
            .map{ taxon, data -> [ taxon, data.value instanceof List ? data.value.taxonomy.classification.species : null ]}
            .filter{ taxon, species -> species }
            .set{ ch_ncbi_species }

        NCBI_DATA
            .out
            .subtype
            .splitCsv(header: false, quote: '"')
            .map{ formatSubtypeData(it) }
            .set{ ch_ncbi_subtype }

        NCBI_DATA
            .out
            .data_reports
            .transpose()
            .map{ taxon, json -> [ taxon, file(json).getSimpleName(), json ] }
            .splitJson()
            .groupTuple(by: [0,1])
            .map{ [ it[0], formatTaxonData(it) ] }
            .combine(ch_ncbi_species, by: 0)
            .map{ taxon, data, species -> data.species = species.findAll{ it -> data.taxIds.contains(it.id) }.name
                                        data.species = data.species[0]
                                        data }
            .map{ [ it.accession, it ] }
            .join( ch_ncbi_subtype.map{ [ it.accession, it.findAll{ item -> item.key != 'accession' } ] }, by: 0, remainder: true )
            .filter{ accession, main, other -> main }
            .map{ accession, main, other -> main + (other ? other : [ subtype: null, genotype: null, segment: null ]) }
            .set{ ch_ncbi_data }

        ch_ncbi_data
            .filter{ it.segment }
            .map{ [ it.taxon, it.segment ] }
            .groupTuple(by: 0)
            .map{ taxon, segments -> [ taxon, segments.flatten().unique() ] }
            .subscribe{ taxon, segments -> println "${taxon} Segment Options: ${segments}" }
        
        ch_input
            .filter{ it.segmentSynomyms }
            .map{   def segmentOptions = [:]
                    it.segmentSynomyms.split(';').each{ s -> def syns = s.split('\\|').toList()
                                                             segmentOptions[ syns.get(0) ] = (syns + syns.collect{ v -> v.toUpperCase() } + syns.collect{ v -> v.toLowerCase() }).unique()
                                                    }
                    [ it.taxon, segmentOptions ]
             }
             .set{ ch_segmentOptions }
        ch_ncbi_data
            .map{ [ it.taxon, it ] }
            .groupTuple(by: 0)
            .join(ch_segmentOptions, by: 0, remainder: true)
            .transpose()
            .map{ fixSegmentSynonyms(it) }
            .groupTuple(by: 0)
            .map{ taxon, data -> [ data,  data.segment.size() > 0 ? true : false ] }
            .transpose()
            .map{ data, status -> data + [ segmented: true ] }
            .filter{ (it.segmented && it.segment) || (! it.segmented) }
            .set{ ch_ncbi_data }
    }

    // /*
    // =============================================================================================================================
    //     QUALITY FILTER SEQUENCES
    // =============================================================================================================================
    // */
    
    // // MODULE: Filter low quality sequences & remove duplicates
    // INPUT_QC(
    //     ch_input.map{ [ it.taxon, it.segment, it.assembly, it.length ] }
    // )
    // ch_versions = ch_versions.mix(INPUT_QC.out.versions.first())


    // /*
    // =============================================================================================================================
    //     CLUSTER SEQUENCES: CLUSTER SUBSET
    // =============================================================================================================================
    // */
    // // MODULE: Determine pairwise Mash distances of sequence subset (number determined by `--max_clusters`; will include all if lower than the assigned threshold.)
    // MASH_TOP (
    //     INPUT_QC.out.seqs.map{taxon, segment, top, remainder, remainder_count -> [ taxon, segment, top ]}
    // )
    // ch_versions = ch_versions.mix(MASH_TOP.out.versions.first())

    // // MODULE: Cluster subset with hclust & cutree
    // CLUSTER_MAIN (
    //     MASH_TOP.out.dist,
    //     "main"
    // )
    // ch_versions = ch_versions.mix(CLUSTER_MAIN.out.versions.first())

    // // Set main cluster channel
    // CLUSTER_MAIN
    //     .out
    //     .results
    //     .map{ taxon, segment, results -> results }
    //     .splitCsv(header: true)
    //     .set{ ch_cluster_main }

    // /*
    // =============================================================================================================================
    //     CLUSTER SEQUENCES: ASSIGN REMAINING USING SUBSET CLASSIFIER
    // =============================================================================================================================
    // */
    // // MODULE: Run Mash on the remainder of sequences compared to representatives of each cluster identified in round 1
    // MASH_REMAINDER (
    //     INPUT_QC
    //         .out
    //         .seqs
    //         .filter{ taxon, segment, top, remainder, remainder_count -> remainder_count.toInteger() > 0 }
    //         .map{ taxon, segment, top, remainder, remainder_count -> [ taxon, segment, top, remainder ] }
    //         .join(CLUSTER_MAIN.out.results, by: [0,1])
    // )
    // // MODULE: Assign the remainder of sequences to a cluster
    // ASSIGN_REMAINDER (
    //     MASH_REMAINDER.out.results
    // )
    // ch_versions = ch_versions.mix(ASSIGN_REMAINDER.out.versions.first())

    // // Set assigned cluster channel
    // ASSIGN_REMAINDER
    //     .out
    //     .assigned
    //     .map{ taxon, segment, results -> results }
    //     .splitCsv( header: true )
    //     .set{ ch_cluster_assigned }

    // /*
    // =============================================================================================================================
    //     CLUSTER SEQUENCES: CLUSTER LOOSE-ENDS
    // =============================================================================================================================
    // */
    // SEQTK_LOOSEENDS (
    //     ASSIGN_REMAINDER
    //         .out
    //         .not_assigned
    //         .filter{ taxon, segment, seq_list, count -> count.toInteger() > 1 }
    //         .map{ taxon, segment, seq_list, count -> [ taxon, segment, seq_list ] }
    //         .join(INPUT_QC.out.seqs.map{ taxon, segment, top, remainder, remainder_count -> [ taxon, segment, remainder ] }, by: [0,1])
    // )
    // ch_versions = ch_versions.mix(SEQTK_LOOSEENDS.out.versions.first())

    // // MODULE: Determine pairwise Mash distances of loose-ends.
    // MASH_LOOSEENDS (
    //     SEQTK_LOOSEENDS.out.sequences
    // )
    // ch_versions = ch_versions.mix(MASH_LOOSEENDS.out.versions.first())

    // // MODULE: Cluster loose-ends with hclust & cutree
    // CLUSTER_LOOSEENDS (
    //     MASH_LOOSEENDS.out.dist,
    //     "looseends"
    // )
    // ch_versions = ch_versions.mix(CLUSTER_LOOSEENDS.out.versions.first())

    // // Set loose-ends cluster channel
    // CLUSTER_LOOSEENDS
    //     .out
    //     .results
    //     .map{ taxon, segment, results -> results }
    //     .splitCsv( header: true )
    //     .set{ ch_cluster_looseends }

    // /*
    // =============================================================================================================================
    //     CLUSTER SEQUENCES: COMBINE ALL CLUSTERS
    // =============================================================================================================================
    // */
    // // Combine main cluster results and assigned cluster results
    // ch_cluster_main
    //     .concat(ch_cluster_assigned)
    //     .set{ ch_clusters }
    // // Adjust loose-end cluster numbers based on the main clustering numbers
    // ch_clusters
    //     .map{ [ it.taxon, it.segment, it ] }
    //     .groupTuple(by: [0,1])
    //     .map{ taxon, segment, data -> [ taxon, segment, data.max{ it.cluster }.cluster ] }
    //     .join( ch_cluster_looseends.map{ [ it.taxon, it.segment, it ] }.groupTuple(by: [0,1]), by: [0,1] )
    //     .transpose()
    //     .map{ taxon, segment, max_cluster, it -> it.cluster = it.cluster.toInteger() + max_cluster.toInteger()
    //                                             it }
    //     .concat(ch_clusters)
    //     .map{ it.cluster = it.cluster.toString()
    //           it }
    //     .set{ ch_clusters }

    // // Save to file
    // ch_clusters
    //     .take(1)
    //     .map{ it.keySet().join(',') }
    //     .concat(ch_clusters.map{ it.values().join(',').replace(' ','_') })
    //     .collectFile(name: 'clusters.csv', newLine: true, sort: 'index')
    //     .set{ ch_clusters_file }

    // // MODULE: Split clusters into multi-fasta files
    // SEQTK_SUBSEQ(
    //     ch_clusters
    //         .map{ [ it.taxon, it.segment, it.cluster, it.seq ] }
    //         .groupTuple(by: [0,1,2])
    //         .combine(INPUT_QC.out.all, by: [0,1])        
    // )
    // ch_versions = ch_versions.mix(SEQTK_SUBSEQ.out.versions.first())

    // /*
    // =============================================================================================================================
    //     ALIGN SEQUENCE CLUSTERS 
    // =============================================================================================================================
    // */
    // // MODULE: Align clustered sequences with mafft - only performed on clusters containing more than one sequence 
    // MAFFT(
    //     SEQTK_SUBSEQ
    //         .out
    //         .sequences
    //         .filter{ taxon, segment, cluster, seqs, n_seq -> n_seq.toInteger() > 1 }
    //         .map{ taxon, segment, cluster, seqs, n_seq -> [ taxon, segment, cluster, seqs ] }
    // )
    // ch_versions = ch_versions.mix(MAFFT.out.versions.first())

    // // recombine with singletons (i.e., clusters containing 1 sequence)
    // SEQTK_SUBSEQ
    //     .out
    //     .sequences
    //     .filter{ taxon, segment, cluster, seqs, n_seq -> n_seq.toInteger() == 1 }
    //     .map{ taxon, segment, cluster, seqs, n_seq -> [ taxon, segment, cluster, seqs ] }
    //     .concat(MAFFT.out.fa)
    //     .set{ ch_alignments }
    // /*
    // =============================================================================================================================
    //     CREATE CONSENSUS 
    // =============================================================================================================================
    // */
    // // MODULE: Create consensus sequences
    // CONSENSUS(
    //     ch_alignments
    // )
    // ch_versions = ch_versions.mix(CONSENSUS.out.versions.first())

    // CONSENSUS.out.fa.groupTuple(by: [0,1]).set{ ch_consensus }

    // /*
    // =============================================================================================================================
    //     CONDENSE CONSENSUS SEQS
    // =============================================================================================================================
    // */
    // // MODULE: Determine pairwise Mash distances of consensus sequences.
    // MASH_CONDENSE (
    //     ch_consensus
    // )
    // ch_versions = ch_versions.mix(MASH_CONDENSE.out.versions.first())

    // // MODULE: Condense sequences that share sequence identity below `--dist_threshold`
    // CONDENSE (
    //     MASH_CONDENSE
    //         .out
    //         .dist
    //         .join(ch_consensus, by: [0,1])
    //         .combine(ch_clusters_file)
    // )
    // ch_versions = ch_versions.mix(CONDENSE.out.versions.first())

    // /*
    // =============================================================================================================================
    //     GATHER DATA ON FINAL SEQUENCES
    // =============================================================================================================================
    // */
    // // MODULE: Determine average nucleotide identity between consensus sequences
    // FASTANI_AVA (
    //     CONDENSE.out.results.map{ taxon, segment, summary, consensus, length -> [ taxon, segment, consensus, length ] }
    // )
    // ch_versions = ch_versions.mix(FASTANI_AVA.out.versions.first())

    // /*
    // =============================================================================================================================
    //     SUMMARIZE RESULTS
    // =============================================================================================================================
    // */
    // // MODULE: Create individual summaries
    // SUMMARY(
    //     CONDENSE
    //         .out
    //         .results
    //         .map{ taxon, segment, summary, assembly, length -> [ taxon, segment, summary ] }
    //         .combine(ch_clusters_file)
    //         .join(INPUT_QC.out.all, by: [0,1])
    //         .join(ch_input.map{ [ it.taxon, it.segment, it.assembly, it.metadata ] }, by: [0,1])
    //         .join(FASTANI_AVA.out.ani, by: [0,1])
    // )
    // ch_versions = ch_versions.mix(SUMMARY.out.versions.first())

    // // MODULE: Export reference samplesheet for VAPER
    // CONDENSE
    //     .out
    //     .results
    //     .map{ taxon, segment, summary, assembly, length -> [ taxon, segment, assembly ] }
    //     .transpose()
    //     .map{ taxon, segment, assembly -> taxon+"\t"+segment+"\t"+assembly.getName() }
    //     .collectFile(name: "sheet.tsv", newLine: true)
    //     .set{ ch_ref_lines }
    // EXPORT (
    //     ch_ref_lines,
    //     ch_timestamp
    // )
    // ch_versions = ch_versions.mix(EXPORT.out.versions.first())
    

    // /*
    // =============================================================================================================================
    //     NEXTFLOW DEFAULTS
    // =============================================================================================================================
    // */
    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
