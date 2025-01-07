#!/usr/bin/env nextflow

params.outdir = params.outdir ? params.outdir : 'results'
if(! params.taxon){ exit 1, "Supply the taxon name using '--taxon' (e.g., Alphainfluenzavirus)" }
if(! params.species){ exit 1, "Supply the species name using '--species' (e.g., Alphainfluenzavirus influenzae)" }
workflow {
    Channel
        .fromPath(params.fasta)
        .splitFasta(record: [id: true, seqString: true])
        .map{ def meta = it.id.split('\\|')
              [ id: meta[0], accession: meta[1], segment_name: meta[3], segment_number: meta[4], seqString: it.seqString ] }
        .set{ ch_fasta }
    ch_fasta
        .map{ [ it.segment_name, it.segment_number, ">${it.accession}\n${it.seqString}" ] }
        .groupTuple(by: [0,1])
        .subscribe{ segment_name, segment_number, records -> def fwork = file(workflow.workDir).resolve("${file(params.fasta).baseName}-seg-${segment_number}-${segment_name}.fa") 
                                                             def fres  = file(params.outdir).resolve("${file(params.fasta).baseName}-seg-${segment_number}-${segment_name}.fa") 
                                                             fwork.text = records.join('\n')
                                                             fwork.copyTo(fres)}

    FORMAT_METADATA(
        file(params.metadata)
    ).splitCsv(header: true)
        .map{ def data = [:]
            it.each{ v -> data[ v.key.replaceAll('"','') ] = v.value ? v.value.replaceAll('"', '') : null } 
            data
        }
        .map{ it.segments = [:]
            it.each{ v -> if( v.key.contains('_Segment_Id') ){
                it.segments[ v.key.replaceAll('_Segment_Id', '') ] = v.value ? v.value.split('\\|')[0] : null
            } }
            def result = [ it.Isolate_Id ] + it.segments + [ 
                                    taxon: params.taxon,
                                    species: params.species,
                                    geographicRegion: it.Location ? it.Location.split(' / ')[0] : null,
                                    collectionDate: it.Collection_Date ? it.Collection_Date.substring(0,4) : null,
                                    organismName_host: it.Host ? it.Host : null,
                                    subtype: it.Subtype ? it.Subtype.replaceAll(/[AB] \/ /, '') : null,
                                    lineage: it.Lineage ? it.Lineage : null ,
                                    clade: it.Clade ? it.Clade : null
                                    ]
        }
        .map{ id, segment_accessions, data -> [ id, segment_accessions, data + splitSubtype( data.subtype ) ]
 }
        .join( ch_fasta.map{ [ it.id, [ it.segment_name, it.segment_number ] ] }.groupTuple( by: 0), by: 0 )
        .map{ id, segment_accessions, data, segment_numbers -> segment_numbers = segment_numbers.collect{ v -> v = v + [ segment_accessions.containsKey(v[0]) ? segment_accessions[v[0]] : null ] }
                                                                [ segment_numbers ] + data  }
        .transpose()
        .map{ segment_info, data -> [ segment_info[0], segment_info[1] ] + [ [ accession: segment_info[2], segment: segment_info[1] ] + data ] }
        .groupTuple(by: [0,1])
        .map{ segment_name, segment_number, data -> def table = channelToTable( data ) 
                                                    def fwork = file(workflow.workDir).resolve("${file(params.fasta).baseName}-seg-${segment_number}-${segment_name}.csv")
                                                    def fres  = file(params.outdir).resolve("${file(params.fasta).baseName}-seg-${segment_number}-${segment_name}.csv")
                                                    fwork.text = table
                                                    fwork.copyTo(fres)
        }
}

def splitSubtype( subtype ){
    def HA_matcher = []
    def NA_matcher = []
    if(subtype){
        HA_matcher = subtype =~ /H[0-9]+/
        NA_matcher = subtype =~ /N[0-9]+/
    }
    
    return [ HA_type: HA_matcher ? HA_matcher[0] : null, NA_type: NA_matcher ? NA_matcher[0] : null ]

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

process FORMAT_METADATA{
    input:
    path metadata

    output:
    path "metadata.formated.csv"

    script:
    """
    Rscript -e 'library(tidyverse); replace_comma <- function(str){ return(str_replace_all(str,",",":")) }; read_csv("${metadata}") %>% select_all(~gsub("\\\\s+|\\\\.", "_", .)) %>% mutate_all(replace_comma) %>% write.csv(file = "metadata.formated.csv", quote = T, row.names = F, na = "")'
    """


}