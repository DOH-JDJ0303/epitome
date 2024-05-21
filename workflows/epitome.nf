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
include { INPUT_QC                     } from '../modules/local/input-qc'
include { MASH_TOP                     } from '../modules/local/mash'
include { CLUSTER                      } from '../modules/local/cluster'
include { MASH_REMAINDER               } from '../modules/local/mash'
include { ASSIGN_REMAINDER             } from '../modules/local/assign-remainder'
include { SEQTK_LOOSEENDS              } from '../modules/local/seqtk_subseq'
include { MASH_TOP as MASH_LOOSEENDS   } from '../modules/local/mash'
include { CLUSTER as CLUSTER_LOOSEENDS } from '../modules/local/cluster'
include { BIND_CLUSTERS                } from '../modules/local/bind-clusters.nf'
include { SEQTK_SUBSEQ                 } from '../modules/local/seqtk_subseq'
include { MAFFT                        } from '../modules/local/mafft'
include { CONSENSUS                    } from '../modules/local/consensus'
include { MASH_TOP as MASH_CONDENSE    } from '../modules/local/mash'
include { CONDENSE                     } from '../modules/local/condense'
include { FASTANI_AVA                  } from '../modules/local/fastani'
include { FASTANI_SEEDS                } from '../modules/local/fastani'
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
        .map{ tuple(it.taxa, it.segment, file(it.assembly, checkIfExists: true), it.length) }
        .set{ manifest }

    /*
    =============================================================================================================================
        QUALITY FILTER SEQUENCES
    =============================================================================================================================
    */
    
    // MODULE: Filter low quality sequences & remove duplicates
    INPUT_QC(
        manifest
    )
    ch_versions = ch_versions.mix(INPUT_QC.out.versions.first())


    /*
    =============================================================================================================================
        CLUSTER SEQUENCES: CLUSTER SUBSET
    =============================================================================================================================
    */
    // MODULE: Determine pairwise Mash distances of sequence subset (number determined by `--max_clusters`; will include all if lower than the assigned threshold.)
    MASH_TOP (
        INPUT_QC.out.seqs.map{taxa, segment, top, remainder, remainder_count -> [ taxa, segment, top ]}
    )
    ch_versions = ch_versions.mix(MASH_TOP.out.versions.first())

    // MODULE: Cluster subset with hclust & cutree
    CLUSTER (
        MASH_TOP.out.dist,
        "main"
    )
    ch_versions = ch_versions.mix(CLUSTER.out.versions.first())

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
            .filter{ taxa, segment, top, remainder, remainder_count -> remainder_count.toInteger() > 0 }
            .map{ taxa, segment, top, remainder, remainder_count -> [ taxa, segment, top, remainder ] }
            .join(CLUSTER.out.results, by: [0,1])
    )
    // MODULE: Assign the remainder of sequences to a cluster
    ASSIGN_REMAINDER (
        MASH_REMAINDER.out.results
    )
    ch_versions = ch_versions.mix(CLUSTER.out.versions.first())

    /*
    =============================================================================================================================
        CLUSTER SEQUENCES: CLUSTER LOOSE-ENDS
    =============================================================================================================================
    */
    SEQTK_LOOSEENDS (
        ASSIGN_REMAINDER
            .out
            .not_assigned
            .filter{ taxa, segment, seq_list, count -> count.toInteger() > 1 }
            .map{ taxa, segment, seq_list, count -> [ taxa, segment, seq_list ] }
            .join(INPUT_QC.out.seqs.map{ taxa, segment, top, remainder, remainder_count -> [ taxa, segment, remainder ] }, by: [0,1])
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

    /*
    =============================================================================================================================
        CLUSTER SEQUENCES: COMBINE ALL CLUSTERS
    =============================================================================================================================
    */
    // MODULE: Combine all cluster channels. Loose-end clusters are adjusted based on largest subset cluster.
    BIND_CLUSTERS (
        CLUSTER.out.results.concat(ASSIGN_REMAINDER.out.assigned).concat(CLUSTER_LOOSEENDS.out.results).groupTuple(by: [0,1])
    )
    ch_versions = ch_versions.mix(BIND_CLUSTERS.out.versions.first())

    BIND_CLUSTERS
        .out
        .results
        .map{ taxa, segment, result -> result }
        .splitCsv(header: true)
        .map{ tuple(it.taxa, it.segment, it.cluster, it.seq) }
        .groupTuple(by: [0,1,2])
        .combine(INPUT_QC.out.all, by: [0,1])
        .map{ taxa, segment, cluster, contigs, seqs -> [ taxa, segment, cluster, contigs, seqs, contigs.size() ] }
        .set{ clusters }
    

    // MODULE: Split clusters into multi-fasta files
    SEQTK_SUBSEQ(
        clusters
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
            .filter{ taxa, segment, cluster, seqs, count -> count > 1 }
            .map{ taxa, segment, cluster, seqs, count -> [ taxa, segment, cluster, seqs ] }
    )
    ch_versions = ch_versions.mix(MAFFT.out.versions.first())

    // recombine with singletons (i.e., clusters containing 1 sequence)
    SEQTK_SUBSEQ
        .out
        .sequences
        .filter{ taxa, segment, cluster, seqs, count -> count == 1 }
        .map{ taxa, segment, cluster, seqs, count -> [ taxa, segment, cluster, seqs ] }
        .concat(MAFFT.out.fa)
        .set{ alignments }
    
    /*
    =============================================================================================================================
        CREATE CONSENSUS 
    =============================================================================================================================
    */
    // MODULE: Create consensus sequences
    CONSENSUS(
        alignments
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
        MASH_CONDENSE.out.dist.join(ch_consensus, by: [0,1]).join(BIND_CLUSTERS.out.results, by: [0,1])
    )
    ch_versions = ch_versions.mix(CONDENSE.out.versions.first())

    /*
    =============================================================================================================================
        GATHER DATA ON FINAL SEQUENCES
    =============================================================================================================================
    */
    // MODULE: Determine average nucleotide identity between consensus sequences
    FASTANI_AVA (
        CONDENSE.out.results.map{ taxa, segment, summary, consensus, length -> [ taxa, segment, consensus, length ] }
    )
    ch_versions = ch_versions.mix(FASTANI_AVA.out.versions.first())
    // Classify consensus sequences based on supplied seed sequences - if supplied
    if(params.seeds){
        Channel
            .fromPath(params.seeds)
            .splitCsv(header:true)
            .map{ tuple(it.ref, file(it.assembly)) }
            .set{ seeds }
        // MODULE: Determine average nucleotide identity between the consensus sequences and seed sequences
        FASTANI_SEEDS (
            CONDENSE.results.map{ taxa, segment, summary, assembly, length -> assembly }.collect(),
            seeds.map{ ref, assembly -> assembly }.collect()
        )    
        ch_versions = ch_versions.mix(FASTANI_SEEDS.out.versions.first())
    }

    /*
    =============================================================================================================================
        SUMMARIZE RESULTS
    =============================================================================================================================
    */
    // MODULE: Create summary
    SUMMARY(
        CONDENSE.out.results.map{ taxa, segment, summary, assembly, length -> summary }.splitText(keepHeader: true).collectFile(name: 'all-summaries.csv'),
        FASTANI_AVA.out.ani.map{ taxa, segment, ani -> ani }.splitText().collectFile(name: "all-ani.tsv"),
        params.seeds ? FASTANI_SEEDS.out.ani : [],
        params.seeds ? file(params.seeds) : [],
        ch_timestamp
    )
    ch_versions = ch_versions.mix(SUMMARY.out.versions.first())

    // MODULE: Export reference samplesheet for VAPER
    CONDENSE
        .out
        .results
        .map{ taxa, segment, summary, assembly, length -> [ taxa, segment, assembly ] }
        .transpose()
        .map{ taxa, segment, assembly -> taxa+"\t"+segment+"\t"+assembly.getName() }
        .collectFile(name: "sheet.tsv", newLine: true)
        .set{ ch_ref_lines }
    EXPORT (
        ch_ref_lines,
        ch_timestamp
    )
    ch_versions = ch_versions.mix(EXPORT.out.versions.first())
    

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
