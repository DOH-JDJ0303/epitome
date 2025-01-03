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
include { INPUT_QC                     } from '../subworkflows/local/input-qc'
include { NCBI_DATA_SUBWF              } from '../subworkflows/local/ncbi-data'

include { TIMESTAMP                    } from '../modules/local/timestamp'
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
        .splitCsv(header:true, quote: '"')
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
        GATHER NCBI DATA
    =============================================================================================================================
    */
    if(params.ncbi){

        NCBI_DATA_SUBWF(
            ch_input
        ).input.set{ ch_input }
    }

    /*
    =============================================================================================================================
        QUALITY FILTER SEQUENCES
    =============================================================================================================================
    */
    INPUT_QC(
        ch_input
    )

    /*
    =============================================================================================================================
        CLUSTER SEQUENCES: CLUSTER SUBSET
    =============================================================================================================================
    */
    // MODULE: Determine pairwise Mash distances of sequence subset (number determined by `--max_clusters`; will include all if lower than the assigned threshold.)
    MASH_TOP (
        INPUT_QC.out.seqs.map{taxon, segment, all, top, remainder, remainder_count -> [ taxon, segment, top ]}
    )
    ch_versions = ch_versions.mix(MASH_TOP.out.versions.first())

    // MODULE: Cluster subset with hclust & cutree
    CLUSTER_MAIN (
        MASH_TOP.out.dist,
        "main"
    )
    ch_versions = ch_versions.mix(CLUSTER_MAIN.out.versions.first())

    // Set main cluster channel
    CLUSTER_MAIN
        .out
        .results
        .map{ taxon, segment, results -> results }
        .splitCsv(header: true, quote: '"')
        .set{ ch_cluster_main }

    /*
    =============================================================================================================================
        CLUSTER SEQUENCES: ASSIGN REMAINING USING SUBSET CLASSIFIER
    =============================================================================================================================
    */
    // MODULE: Run Mash on the remainder of sequences compared to representatives of each cluster identified in round 1
    MASH_REMAINDER (
        INPUT_QC
            .out
            .seqs
            .filter{ taxon, segment, all, top, remainder, remainder_count -> remainder_count.toInteger() > 0 }
            .map{ taxon, segment, all, top, remainder, remainder_count -> [ taxon, segment, top, remainder ] }
            .join(CLUSTER_MAIN.out.results, by: [0,1])
    )
    // MODULE: Assign the remainder of sequences to a cluster
    ASSIGN_REMAINDER (
        MASH_REMAINDER.out.results
    )
    ch_versions = ch_versions.mix(ASSIGN_REMAINDER.out.versions.first())

    // Set assigned cluster channel
    ASSIGN_REMAINDER
        .out
        .assigned
        .map{ taxon, segment, results -> results }
        .splitCsv( header: true, quote: '"' )
        .set{ ch_cluster_assigned }

    /*
    =============================================================================================================================
        CLUSTER SEQUENCES: CLUSTER LOOSE-ENDS
    =============================================================================================================================
    */
    SEQTK_LOOSEENDS (
        ASSIGN_REMAINDER
            .out
            .not_assigned
            .filter{ taxon, segment, seq_list, count -> count.toInteger() > 1 }
            .map{ taxon, segment, seq_list, count -> [ taxon, segment, seq_list ] }
            .join(INPUT_QC.out.seqs.map{ taxon, segment, all, top, remainder, remainder_count -> [ taxon, segment, remainder ] }, by: [0,1])
    )
    ch_versions = ch_versions.mix(SEQTK_LOOSEENDS.out.versions.first())

    // MODULE: Determine pairwise Mash distances of loose-ends.
    MASH_LOOSEENDS (
        SEQTK_LOOSEENDS.out.sequences
    )
    ch_versions = ch_versions.mix(MASH_LOOSEENDS.out.versions.first())

    // MODULE: Cluster loose-ends with hclust & cutree
    CLUSTER_LOOSEENDS (
        MASH_LOOSEENDS.out.dist,
        "looseends"
    )
    ch_versions = ch_versions.mix(CLUSTER_LOOSEENDS.out.versions.first())

    // Set loose-ends cluster channel
    CLUSTER_LOOSEENDS
        .out
        .results
        .map{ taxon, segment, results -> results }
        .splitCsv( header: true, quote: '"' )
        .set{ ch_cluster_looseends }

    /*
    =============================================================================================================================
        CLUSTER SEQUENCES: COMBINE ALL CLUSTERS
    =============================================================================================================================
    */
    // Combine main cluster results and assigned cluster results
    ch_cluster_main
        .concat(ch_cluster_assigned)
        .set{ ch_clusters }
    // Adjust loose-end cluster numbers based on the main clustering numbers
    ch_clusters
        .map{ [ it.taxon, it.segment, it ] }
        .groupTuple(by: [0,1])
        .map{ taxon, segment, data -> [ taxon, segment, data.max{ it.cluster }.cluster ] }
        .join( ch_cluster_looseends.map{ [ it.taxon, it.segment, it ] }.groupTuple(by: [0,1]), by: [0,1] )
        .transpose()
        .map{ taxon, segment, max_cluster, it -> it.cluster = it.cluster.toInteger() + max_cluster.toInteger()
                                                it }
        .concat(ch_clusters)
        .map{ it.cluster = it.cluster.toString()
              it }
        .set{ ch_clusters }

    // Save to file
    ch_clusters
        .map{ [ it.taxon, it.segment, it ] }
        .groupTuple(by: [0,1])
        .map{ taxon, segment, data -> def table = channelToTable(data)
                                      def fwork = file(workflow.workDir).resolve("${taxon}-${segment}-clusters.csv")
                                      def fres  = file(params.outdir).resolve(taxon).resolve(segment).resolve('clusters').resolve("clusters.all.csv")
                                      fwork.text = table
                                      fwork.copyTo(fres)
                                      [ taxon: taxon, segment: segment, file: fwork ]
        }
        .set{ ch_clusters_file }
    // MODULE: Split clusters into multi-fasta files
    SEQTK_SUBSEQ(
        ch_clusters
            .map{ [ it.taxon, it.segment, it.cluster, it.seq ] }
            .groupTuple(by: [0,1,2])
            .combine(INPUT_QC.out.seqs.map{ taxon, segment, all, top, remainder, remainder_count -> [ taxon, segment, all ] }, by: [0,1])        
    )
    ch_versions = ch_versions.mix(SEQTK_SUBSEQ.out.versions.first())

    /*
    =============================================================================================================================
        ALIGN SEQUENCE CLUSTERS 
    =============================================================================================================================
    */
    // MODULE: Align clustered sequences with mafft - only performed on clusters containing more than one sequence 
    MAFFT(
        SEQTK_SUBSEQ
            .out
            .sequences
            .filter{ taxon, segment, cluster, seqs, n_seq -> n_seq.toInteger() > 1 }
            .map{ taxon, segment, cluster, seqs, n_seq -> [ taxon, segment, cluster, seqs ] }
    )
    ch_versions = ch_versions.mix(MAFFT.out.versions.first())

    // recombine with singletons (i.e., clusters containing 1 sequence)
    SEQTK_SUBSEQ
        .out
        .sequences
        .filter{ taxon, segment, cluster, seqs, n_seq -> n_seq.toInteger() == 1 }
        .map{ taxon, segment, cluster, seqs, n_seq -> [ taxon, segment, cluster, seqs ] }
        .concat(MAFFT.out.fa)
        .set{ ch_alignments }
    /*
    =============================================================================================================================
        CREATE CONSENSUS 
    =============================================================================================================================
    */
    // MODULE: Create consensus sequences
    CONSENSUS(
        ch_alignments
    )
    ch_versions = ch_versions.mix(CONSENSUS.out.versions.first())

    CONSENSUS.out.fa.groupTuple(by: [0,1]).set{ ch_consensus }

    /*
    =============================================================================================================================
        CONDENSE CONSENSUS SEQS
    =============================================================================================================================
    */
    // MODULE: Determine pairwise Mash distances of consensus sequences.
    MASH_CONDENSE (
        ch_consensus
    )
    ch_versions = ch_versions.mix(MASH_CONDENSE.out.versions.first())

    // MODULE: Condense sequences that share sequence identity below `--dist_threshold`
    CONDENSE (
        MASH_CONDENSE
            .out
            .dist
            .join(ch_consensus, by: [0,1])
            .join(ch_clusters_file.map{ [ it.taxon, it.segment, it.file ] }, by: [0,1])
    )
    ch_versions = ch_versions.mix(CONDENSE.out.versions.first())

    /*
    =============================================================================================================================
        SUMMARIZE RESULTS
    =============================================================================================================================
    */
    // Load metadata & save to file
    ch_input
        .map{ it.metadata }
        .splitCsv(header: true, quote: '"')
        .set{ ch_metadata }
    ch_metadata
        .map{ [ it.taxon, it.segment, it ] }
        .groupTuple(by: [0,1])
        .subscribe{ taxon, segment, data -> def table = channelToTable(data)
                                      def fwork = file(workflow.workDir).resolve("${taxon}-${segment}-metadata.csv")
                                      def fres  = file(params.outdir).resolve(taxon).resolve(segment).resolve('metadata').resolve("metadata.csv")
                                      fwork.text = table
                                      fwork.copyTo(fres)
        }
    // Combine data channels
    INPUT_QC
        .out
        .qc
        .map{ [ it.accessions, it ] }
        .transpose()
        .join(ch_metadata.map{ [ it.accession, it ] }, by: 0, remainder: true) // Input QC + Metadata
        .filter{ accession, qc, meta -> qc }
        .map{ accession, qc, meta -> qc + (meta ? meta : [ metadata: 'none' ]) }
        .map{ [ it.taxon, it.segment, it.seq.toString(), it ] }
        .combine(ch_clusters.map{ [ it.taxon, it.segment, it.seq, it.cluster ] }, by: [0,1,2]) // + Clusters     
        .groupTuple(by: [0,1,4])
        .map{ collapseKeys(it) } // Consolidate rows for groups of taxon, segment, cluster
        .map{ [ it.taxon, it.segment, it.cluster, it ] }
        .join(CONDENSE.out.results.map{ taxon, segment, summary, assembly -> summary }.splitCsv(header: true, quote: '"').map{ [ it.taxon, it.segment, it.cluster, it ] }, by: [0,1,2], remainder: true) // + Final references
        .map{ taxon, segment, cluster, data1, data2 -> def data = data1 + data2
                                                       def unwantedKeys = ['accessions','illegalBases','filters','filter']
                                                       unwantedKeys.each{ key -> data.remove(key) } // Remove unwanted keys
                                                       [ taxon, segment, cluster, data ] } 
        .groupTuple(by: [0,1])
        .subscribe{ taxon, segment, cluster, data -> def table = channelToTable(data)
                                                     def fwork = file(workflow.workDir).resolve("${taxon}-${segment}-summary.csv")
                                                     def fres  = file(params.outdir).resolve(taxon).resolve(segment).resolve("summary.csv")
                                                     fwork.text = table
                                                     fwork.copyTo(fres)
        }
    

    // // MODULE: Export reference samplesheet for VAPER
    CONDENSE
        .out
        .results
        .map{ taxon, segment, summary, assembly -> [ taxa: taxon, segment: segment, assembly: file(params.outdir).resolve(taxon).resolve(segment).resolve('consensus').resolve("${assembly.name}.fa.gz") ] }
        .collect()
        .combine(ch_timestamp)
        .map{ data, timestamp -> def table = channelToTable(data) 
                                 file(workflow.workDir).resolve("${timestamp}-epitome.csv").text = table
        }

    

    /*
    =============================================================================================================================
        NEXTFLOW DEFAULTS
    =============================================================================================================================
    */
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def collapseKeys (row) {
    data = row[3]
    results = [taxon: row[0], segment: row[1], cluster: row[4]]
    data.each{ it.each{ v -> if(v.key != 'seqString' & v.value != 'null' ){ results[v.key] = results.containsKey(v.key) ? ( results[v.key].contains(v.value) ? results[v.key] : results[v.key] + v.value) : [ v.value ] } } }
    return results
}

def channelToTable ( data ){
    // Gather all keys
    def allKeys = []
    data.each{ allKeys = allKeys + it.keySet().toList() }
    allKeys.unique()
    // Create table
    def table = [allKeys.join(',')]
    data.each{ table = table + [ allKeys.collect{ k -> it.containsKey(k) ? "\"${it[k]}\"" : 'null' }.join(',') ] }
    return table.join('\n')
}

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
