//
// Check input samplesheet and get read channels
//

include { NCBI_DATA    } from '../../modules/local/ncbi-data'

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
            .map{ [ it.taxon, it.segmentSynonyms, it.segmented ] }
            .groupTuple(by: 0)
            .map{ taxon, segmentSynonyms, segmented -> [ taxon, segmentSynonyms.first(), segmented.first()  ] }
    )
    ch_versions = ch_versions.mix(NCBI_DATA.out.versions.first())

    NCBI_DATA
        .out
        .man
        .splitCsv(header: true)
        .map{taxon, csv -> [csv.prefix, taxon, csv.segment]}
        .set{ ch_ncbi_man }

    NCBI_DATA
        .out
        .fa
        .transpose()
        .map{ taxon, f -> [file(f).getName().replace('.fa.gz', ''), f] }
        .join(ch_ncbi_man, by: 0)
        .map{ prefix, f, taxon, segment -> [taxon, segment, f] }
        .set{ ch_ncbi_fa }

    NCBI_DATA
        .out
        .json
        .transpose()
        .map{ taxon, f -> [file(f).getName().replace('.json', ''), f] }
        .join(ch_ncbi_man, by: 0)
        .map{ prefix, f, taxon, segment -> [taxon, segment, f] }
        .set{ ch_ncbi_json }

    ch_ncbi_fa
        .join(ch_ncbi_json, by: [0,1])
        .map{ taxon, segment, assembly, metadata -> [taxon: taxon, segment: segment, assembly: assembly, metadata: metadata] }
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