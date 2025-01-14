//
// Check input samplesheet and get read channels
//

include { INPUT_QC                     } from '../../modules/local/input-qc'
include { MASH_TOP                     } from '../../modules/local/mash'
include { CLUSTER as CLUSTER_MAIN      } from '../../modules/local/cluster'
include { MASH_REMAINDER               } from '../../modules/local/mash'
include { ASSIGN_REMAINDER             } from '../../modules/local/assign-remainder'
include { SEQTK_LOOSEENDS              } from '../../modules/local/seqtk_subseq'
include { MASH_TOP as MASH_LOOSEENDS   } from '../../modules/local/mash'
include { CLUSTER as CLUSTER_LOOSEENDS } from '../../modules/local/cluster'
include { MERGE_CLUSTERS               } from '../../modules/local/merge-clusters'
include { SEQTK_SUBSEQ                 } from '../../modules/local/seqtk_subseq'
include { MAFFT                        } from '../../modules/local/mafft'
include { CONSENSUS                    } from '../../modules/local/consensus'
include { MASH_TOP as MASH_CONDENSE    } from '../../modules/local/mash'
include { CONDENSE                     } from '../../modules/local/condense'
include { SUMMARY                      } from '../../modules/local/summary'

workflow CREATE_SUBWF {
    take:
    ch_input // [ va(taxon), val(segment), path(assembly), ... ]

    main:
    ch_versions = Channel.empty()
    /*
    =============================================================================================================================
        QUALITY FILTER SEQUENCES
    =============================================================================================================================
    */
    INPUT_QC(
        ch_input
            .map{ [ it.taxon, it.segment, it.assembly, it.exclusions ] }
    )
    ch_versions = ch_versions.mix(INPUT_QC.out.versions.first())

    INPUT_QC
        .out
        .seqs
        .map{taxon, segment, all, top, remainder, status -> [ taxon: taxon, segment: segment, all: all, top: status == 'null' ? null : top, remainder: status == 'null' ? null : remainder, status: status == 'null' ? null : status ]}
        .set{ ch_input_qc }

    /*
    =============================================================================================================================
        CLUSTER SEQUENCES: CLUSTER SUBSET
    =============================================================================================================================
    */
    // MODULE: Determine pairwise Mash distances of sequence subset (number determined by `--max_clusters`; will include all if lower than the assigned threshold.)
    MASH_TOP (
        ch_input_qc
            .map{ [ it.taxon, it.segment, it.status ? it.top : it.all ] }
    )
    ch_versions = ch_versions.mix(MASH_TOP.out.versions.first())

    // MODULE: Cluster subset with hclust & cutree
    CLUSTER_MAIN (
        MASH_TOP.out.dist,
        "main"
    )
    ch_versions = ch_versions.mix(CLUSTER_MAIN.out.versions.first())

    /*
    =============================================================================================================================
        CLUSTER SEQUENCES: ASSIGN REMAINING USING SUBSET CLASSIFIER
    =============================================================================================================================
    */
    // MODULE: Run Mash on the remainder of sequences compared to representatives of each cluster identified in round 1
    MASH_REMAINDER (
        ch_input_qc
            .filter{ it.status }
            .map{ [ it.taxon, it.segment, it.top, it.remainder ] }
            .join(CLUSTER_MAIN.out.results, by: [0,1])
    )
    // MODULE: Assign the remainder of sequences to a cluster
    ASSIGN_REMAINDER (
        MASH_REMAINDER.out.results
    )
    ch_versions = ch_versions.mix(ASSIGN_REMAINDER.out.versions.first())

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
            .join(ch_input_qc.map{ [ it.taxon, it.segment, it.remainder ] }, by: [0,1])
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
    MERGE_CLUSTERS (
        CLUSTER_MAIN
            .out
            .results
            .join(ASSIGN_REMAINDER.out.assigned, by: [0,1], remainder: true)
            .join(CLUSTER_LOOSEENDS.out.results, by: [0,1], remainder: true)
            .map{ taxon, segment, top, assigned, remainder -> [ taxon, segment, top, assigned ? assigned : [], remainder ? remainder : [] ] }
    )
    MERGE_CLUSTERS
        .out
        .merged
        .map{ taxon, cluster, file -> file }
        .splitCsv(header: true, quote: '"')
        .set{ ch_clusters }
    // MODULE: Split clusters into multi-fasta files
    SEQTK_SUBSEQ(
        ch_clusters
            .map{ [ it.taxon, it.segment, it.cluster, it.seq ] }
            .groupTuple(by: [0,1,2])
            .combine(ch_input_qc.map{ [ it.taxon, it.segment, it.all ] }, by: [0,1])        
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
            .join(MERGE_CLUSTERS.out.merged, by: [0,1])
    )
    ch_versions = ch_versions.mix(CONDENSE.out.versions.first())

    /*
     =============================================================================================================================
         SUMMARIZE RESULTS
     =============================================================================================================================
    */
    SUMMARY (
        INPUT_QC
            .out
            .summary
            .join(MERGE_CLUSTERS.out.merged, by: [0,1])
            .join(CONDENSE.out.results.map{ taxon, segment, summary, assembly -> [ taxon, segment, summary ] }, by: [0,1])
            .join( ch_input.map{ [ it.taxon, it.segment, it.metadata ] }, by: [0,1] )
    )

    emit:
    versions = ch_versions
}

/*
=============================================================================================================================
    FUNCTIONS
=============================================================================================================================
*/
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