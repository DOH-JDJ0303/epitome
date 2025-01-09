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
    /*
    =============================================================================================================================
        GATHER DATA FROM NCBI
    =============================================================================================================================
    */
    // MODULE: Gather NCBI data for taxon
    NCBI_DATA(
        ch_input.map{ it.taxon }.unique()
    )
    // Load detailed taxonomy data and extract species value
    NCBI_DATA
        .out
        .taxids
        .splitJson()
        .map{ taxon, data -> [ taxon, data.value instanceof List ? data.value.taxonomy.classification.species : null ]}
        .filter{ taxon, species -> species }
        .set{ ch_ncbi_species }
    // Load NCBI subtype data and extract values (subtype, genotype, and segment)
    NCBI_DATA
        .out
        .subtype
        .splitCsv(header: false, quote: '"')
        .map{ formatSubtypeData(it) }
        .set{ ch_ncbi_subtype }
    // Load limited sample data from NCBI datasets
    NCBI_DATA
        .out
        .data_reports
        .transpose()
        .map{ taxon, json -> [ taxon, file(json).getSimpleName(), json ] }
        .splitJson()
        .groupTuple(by: [0,1])
        .map{ [ it[0], formatTaxonData(it) ] }
        .set{ ch_ncbi_datasets }
    // Combine datastreams
    ch_ncbi_datasets
        .combine(ch_ncbi_species, by: 0)
        .map{ taxon, data, species -> data.species = species.findAll{ it -> data.taxIds.contains(it.id) }.name
                                    data.species = data.species[0]
                                    data }
        .map{ [ it.accession, it ] }
        .join( ch_ncbi_subtype.map{ [ it.accession, it.findAll{ item -> item.key != 'accession' } ] }, by: 0, remainder: true )
        .filter{ accession, main, other -> main }
        .map{ accession, main, other -> main + (other ? other : [ subtype: null, genotype: null, segment: null ]) }
        .set{ ch_ncbi_data }
    /*
    =============================================================================================================================
        CLEAN UP SEGMENTS
    =============================================================================================================================
    */
    // Report the segment options per taxon (for troubleshooting purposes)
    ch_ncbi_data
        .filter{ it.segment }
        .map{ [ it.taxon, it.segment ] }
        .groupTuple(by: 0)
        .map{ taxon, segments -> [ taxon, segments.flatten().unique() ] }
        .subscribe{ taxon, segments -> println "${taxon} Segment Options: ${segments}" }
    // Parse segment synonyms (if supplied)
    ch_input
        .filter{ it.segmentSynomyms }
        .map{   def segmentOptions = [:]
                it.segmentSynomyms.split(';').each{ s -> def syns = s.split('\\|').toList()
                                                         segmentOptions[ syns.get(0) ] = (syns + syns.collect{ v -> v.toUpperCase() } + syns.collect{ v -> v.toLowerCase() }).unique()
                                                }
                [ it.taxon, segmentOptions ]
         }
        .unique()
        .set{ ch_segmentOptions }
    // Standardize segment names based on supplied synonyms
    ch_ncbi_data
        .map{ [ it.taxon, it ] }
        .groupTuple(by: 0)
        .join(ch_segmentOptions, by: 0, remainder: true)
        .transpose()
        .map{ fixSegmentSynonyms(it) }
        .groupTuple(by: 0)
        .map{ taxon, data -> [ data,  data.findAll{ it.segment } ? true : false ] }
        .transpose()
        .map{ data, status -> data + [ segmented: status ] }
        .filter{ (it.segmented && it.segment) || (! it.segmented) }
        .map{ it.segment = it.segmented ? it.segment : 'wg'
              it.remove('taxIds')
              it }
        .set{ ch_ncbi_data }
    ch_ncbi_data
        .map{ [ it.taxon, it.segment, it ] }
        .groupTuple(by: [0,1])
        .map{ taxon, segment, data -> def table = channelToTable(data)
                                      def fwork = file(workflow.workDir).resolve("${taxon}-${segment}-ncbi-data.csv")
                                      fwork.text = table
                                      [ taxon: taxon, segment: segment, file: fwork ]
        }
        .set{ch_ncbi_data_file}
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
            .map{ it } 
            .concat(ch_input.filter{ it.segment && it.assembly }.map{ it })
            .map{ [ it.taxon, it.segment, it.assembly, it.metadata ] }
            .groupTuple(by: [0,1])
            .map{ taxon, segment, assembly, metadata -> [ taxon, segment, assembly, metadata ] }
    )

    emit:
    input = MERGE_INPUTS.out.merged.map{ taxon, segment, assembly, metadata -> [ taxon: taxon, segment: segment, assembly: assembly, metadata: metadata ] }
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

def formatTaxonData (row) {
    def data = [:]
    row[2].each{ data[it.key] = it.value }

    def result = [ taxon: row[0], 
                   accession: data.accession, 
                   collectionDate: data.containsKey('isolate') ? ( data.isolate.containsKey('collectionDate') ? data.isolate.collectionDate : null ) : null, 
                   geographicRegion: data.containsKey('location') ? ( data.location.containsKey('geographicRegion') ? data.location.geographicRegion : null ) : null,
                   organismName_host: data.containsKey('host') ? ( data.host.containsKey('organismName') ? data.host.organismName : null ) : null,
                   organismName_virus: data.virus.organismName.split(' ').findAll{ ! it.contains('/') }.join(' '),
                   taxIds: data.virus.lineage.taxId
                   ]
    return result
}

def formatSubtypeData (row) {
    def target_keys = ['segment','subtype','genotype','serotype']
    def results     = [ accession: row[1][0] ]
    def keys        = row[1][1].split('\\|').toList()
    def values      = row[1][2].split('\\|').toList()
    if( keys.size() == values.size() ){
        values.eachWithIndex{ value, index -> results[ keys[index] ] = value }
    }
    target_keys.findAll{ ! results.keySet().contains(it) }.each{ results[it] = null }

    return  results.findAll { it -> ( ['accession'] + target_keys ).contains(it.key) }
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