//
// Check input samplesheet and get read channels
//

include { INPUT_QC       } from '../../modules/local/input-qc'

workflow INPUT_QC_SUBWF {
    take:
    ch_input // [ va(taxon), val(segment), path(assembly), ... ]

    main:


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