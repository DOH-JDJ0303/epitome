//
// Check input samplesheet and get read channels
//

include { INPUT_QC        } from '../../modules/local/input-qc'
include { CLUSTER         } from '../../modules/local/cluster'
include { SEQTK_LOOSEENDS } from '../../modules/local/seqtk_subseq'
include { MAFFT           } from '../../modules/local/mafft'
include { CONSENSUS       } from '../../modules/local/consensus'
include { CONDENSE        } from '../../modules/local/condense'
include { SUMMARY         } from '../../modules/local/summary'

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
    INPUT_QC.out.seqs.set{ ch_input_qc }

    CLUSTER (
        ch_input_qc
    )
    CLUSTER.out.json.set{ ch_clusters }

    /*
    =============================================================================================================================
        ALIGN SEQUENCE CLUSTERS 
    =============================================================================================================================
    */
    // MODULE: Align clustered sequences with mafft - only performed on clusters containing more than one sequence 
    MAFFT(
        CLUSTER
            .out
            .multi
            .transpose()
            .map{ taxon, segment, seqs -> [ taxon, segment, seqs.simpleName.tokenize('-').last(), seqs ] }
    )
    ch_versions = ch_versions.mix(MAFFT.out.versions.first())

    // recombine with singletons (i.e., clusters containing 1 sequence)
    CLUSTER
        .out
        .single
        .transpose()
        .map{ taxon, segment, seqs -> [ taxon, segment, seqs.simpleName.tokenize('-').last(), seqs ] }
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
    // MODULE: Condense sequences that share sequence identity below `--dist_threshold`
    CONDENSE (
        ch_consensus.join(ch_clusters, by: [0,1])
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
            .join(ch_clusters, by: [0,1])
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