//
// Check input samplesheet and get read channels
//

include { NCBI_DATA    } from '../../modules/local/ncbi-data'
include { SEQTK_NCBI   } from '../../modules/local/seqtk_subseq'
include { MERGE_INPUTS } from '../../modules/local/merge-inputs'

workflow NCBI_DATA_SUBWF {
    take:
    ch_input // [ va(taxon), val(segment), path(assembly), ... ]

    main:
    ch_versions = Channel.empty()

    /*
    =============================================================================================================================
        GATHER DATA FROM NCBI
    =============================================================================================================================
    */
    // MODULE: Gather NCBI data for taxon
    NCBI_DATA(
        ch_input
            .map{ [ it.taxon, it.segmentSynomyms, it.segmented ] }
            .groupTuple(by: 0)
            .map{ taxon, segmentSynonyms, segmented -> [ taxon, segmentSynonyms.first(), segmented.first()  ] }
    )
    NCBI_DATA
      .out
      .complete
      .splitCsv(header: true, quote: '"')
      .map{ taxon, data -> data + [ taxon: taxon ]}
      .set{ ch_ncbi_data }
    ch_ncbi_data.filter{it.segment == null}.view()
    ch_ncbi_data
        .map{ [ it.taxon, it.segment, it ] }
        .groupTuple(by: [0,1])
        .map{ taxon, segment, data -> def table = channelToTable(data)
                                      def fwork = file(workflow.workDir).resolve("${taxon}-${segment}-ncbi-data.csv")
                                      fwork.text = table
                                      [ taxon: taxon, segment: segment, file: fwork ] }
        .set{ ch_ncbi_data_file }
    

    /*
    =============================================================================================================================
        PARSE SEQUENCES BY TAXON AND SEGMENT
    =============================================================================================================================
    */
    // Combine data with sequences
    SEQTK_NCBI(
        NCBI_DATA
            .out
            .genomic
            .combine(ch_ncbi_data.map{ [ it.taxon, it.segment, it.accession ] }.groupTuple(by: [0,1] ), by: 0)
    )
    SEQTK_NCBI
        .out
        .sequences
        .join(ch_ncbi_data_file.map{ [ it.taxon, it.segment, it.file ] }, by: [0,1])
        .map{ taxon, segment, assembly, metadata -> [ taxon: taxon, segment: segment, assembly: assembly, metadata: metadata ] }
        .set{ch_input_ncbi}
    /*
    =============================================================================================================================
        MERGE METADATA
    =============================================================================================================================
    */
    // Merge inputs (from NCBI and from the user)
    MERGE_INPUTS (
        ch_input_ncbi
            .concat(ch_input.filter{ it.segment && it.assembly })
            .map{ [ it.taxon, it.segment, it.assembly, it.metadata ] }
            .groupTuple(by: [0,1])
            .map{ taxon, segment, assembly, metadata -> [ taxon, segment, assembly, metadata ] }
    )

    MERGE_INPUTS
        .out
        .merged
        .join( ch_input.map{ [ it.taxon, it.segment, it.exclusions ] }, by: [0,1], remainder: true)
        .map{ taxon, segment, assembly, metadata, exclusions -> [ taxon: taxon, segment: segment, assembly: assembly, metadata: metadata, exclusions: exclusions ? exclusions : [] ] }
        .set{ ch_input }

    emit:
    input    = ch_input
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

def fixSegmentSynonyms(row){
    def taxon = row[0]
    def data  = row[1]
    def syns  = row[2]
    if(syns){
        data.segment = syns.find { key, value -> value.any { it == data.segment }}?.key
    }

    return [ taxon, data ]
}