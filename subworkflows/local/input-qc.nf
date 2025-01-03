//
// Check input samplesheet and get read channels
//

include { EXPORT_QC_SEQS } from '../../modules/local/export-qc-seqs'    

workflow INPUT_QC {
    take:
    ch_input // [ va(taxon), val(segment), path(assembly), ... ]

    main:

    // Gather QC metrics
    ch_input
        .map{ [ [ taxon: it.taxon, segment: it.segment ], it.assembly ] }
        .splitFasta(elem: 1, record: [id: true, seqString: true])
        .groupTuple(by: 0)
        .map{ gatherQcMetrics( it, params.len_threshold, params.amb_threshold ) }
        .flatten()
        .set{ch_qc}

    // Save QC metrics to file
    ch_qc
        .map{ [ it.taxon, it.segment, it ] }
        .groupTuple(by: [0,1])
        .map{ taxon, segment, data -> [ taxon, segment, channelToTable(data) ] }
        .subscribe{ taxon, segment, table -> def f = file(workflow.workDir).resolve("${taxon}-${segment}-input-qc.csv")
                                             f.text = table
                                             f.copyTo(file(params.outdir).resolve(taxon).resolve(segment).resolve('qc').resolve('input-qc.csv'))
        }

    // Export passing sequences to FASTQ
    EXPORT_QC_SEQS(
        ch_qc
            .filter{ ! it.filter }
            .map{ [ it.taxon, it.segment, "${it.seq}\t${it.seqString}" ] }
            .groupTuple(by: [0,1])
            .map{ taxon, segment, records -> file(workflow.workDir).resolve("${taxon}-${segment}-records.txt").text = records.join('\n')
                                             [ taxon, segment, file(workflow.workDir).resolve("${taxon}-${segment}-records.txt") ] }
    )

    emit:
    qc   = ch_qc
    seqs = EXPORT_QC_SEQS.out.seqs
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


def gatherQcMetrics (row, len_threshold, amb_threshold){
    // Accepted bases
    def legalBases = '-ATCGRYSWKMBDHVN'
    // Parse inputs
    meta = row[0]
    records = row[1]
    // Collapse duplicates and gather metrics
    def result = records
                    .groupBy{ it.seqString }
                    .collect{ key, value -> [ taxon: meta.taxon,
                                              segment: meta.segment,
                                              seqString: key, 
                                              accessions: value.id,
                                              length: key.length(),
                                              illegalBases: key.replaceAll( legalBases, '' ).split('.').toList().unique(),
                                              ambRatio: (key =~ 'N-').size() / key.length(),
                                              filters: [:]
                                             ] 
    }
    // Assign sequence numbers
    result.eachWithIndex{ value, index -> result[index]['seq'] = index + 1 }

    // Filter 1: Sequence contains illegal bases
    result = result.collect{ it.filters['illegalBases'] = it.illegalBases ? 'fail' : 'pass'
                             it 
    }
    // Filter 2: Sequence contains too many ambiguous bases
    result = result.collect{ it.filters['ambRatio'] = it.ambRatio >= amb_threshold ? 'fail' : 'pass'
                             it 
    }
    // Filter 3: Sequence length is +/- ratio of median
    def length_med = result.findAll{ it -> it.seq == (( result.seq.max() - 1) / 2).round() }.length.toList().get(0)
    result = result.collect{ it.filters['length'] = it['length'] >= (1 + len_threshold)*length_med || it['length'] <= (1 - len_threshold)*length_med ? 'fail' : 'pass'
                             it 
    }

    // Apply all filters
    result = result.collect { it.filter = it.filters.values().contains('fail') 
                              it
    }

    return result
}
