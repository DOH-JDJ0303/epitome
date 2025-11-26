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
include { NCBI_DATA_SUBWF } from '../subworkflows/local/ncbi-data'
include { MERGE_INPUTS    } from '../modules/local/merge-inputs'
include { CREATE_SUBWF    } from '../subworkflows/local/create'

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
        LOAD SAMPLESHEET
    =============================================================================================================================
    */
    // Load samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(header:true, quote: '"')
        .map{ [ taxon:    it.taxon, 
                assembly: it.containsKey('assembly') ? ( it.assembly ? it.assembly.split(';').collect{f -> file(f, checkIfExists: true)} : [] ) : [], 
                metadata: it.containsKey('metadata') ? ( it.metadata ? it.metadata.split(';').collect{f -> file(f, checkIfExists: true)} : [] ) : [],
                segmented: it.containsKey('segmented') ? it.segmented.toLowerCase() == 'true' : false,
                exclusions: it.containsKey('exclusions') ? ( it.exclusions ? file(it.exclusions, checkIfExists: true) : [] ) : [],
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
        )
        ch_versions = ch_versions.mix(NCBI_DATA_SUBWF.out.versions)
        
        ch_input
            .concat( NCBI_DATA_SUBWF.out.input )
            .map{ [it.taxon, it] }
            .groupTuple()
            .map { String taxon, List<Map> recs ->

                // Use insertion-ordered sets for unique atomic values
                def acc = [:].withDefault { new LinkedHashSet() }

                recs.each { rec ->
                    rec.each { k, v ->
                        if (v == null) return

                        // --- FLATTEN ANY NESTED LISTS / TUPLES ---
                        def flatVals
                        if (v instanceof Collection || v instanceof Object[] ) {
                            flatVals = v.flatten().findAll { it != null }
                        } else {
                            flatVals = [v]
                        }

                        flatVals.each { acc[k] << it }
                    }
                }

                // collapse singletons; multi-value stays a simple list
                def collapsed = acc.collectEntries { k, set ->
                    def vals = set as List
                    [(k): (vals.size() == 1 ? vals[0] : vals)]
                }

                // explicitly retain grouping key
                [taxon: taxon] + collapsed
            }
            .set { ch_input }

    }

    /*
    =============================================================================================================================
        MERGE INPUTS
    =============================================================================================================================
    */
    // Reformat and merge inputs (from NCBI and the user)
    MERGE_INPUTS (
        ch_input.map{ [ it.taxon, it.assembly, it.metadata, it.segmented ] }
    )
    ch_versions = ch_versions.mix(MERGE_INPUTS.out.versions.first())
    MERGE_INPUTS
        .out
        .man
        .splitCsv(header: true)
        .set{ch_input_man}
    
    MERGE_INPUTS
        .out
        .fa
        .transpose()
        .map{ taxon, f -> [taxon, file(f).getName(), f] }
        .join(ch_input_man.map{ taxon, data -> [taxon, data.assembly, data.segment] }, by: [0,1])
        .map{ taxon, f_name, f, segment -> [taxon, segment, f] }
        .set{ch_input_fa}

    MERGE_INPUTS
        .out
        .meta
        .transpose()
        .map{ taxon, f -> [taxon, file(f).getName(), f] }
        .join(ch_input_man.map{ taxon, data -> [taxon, data.metadata, data.segment] }, by: [0,1])
        .map{ taxon, f_name, f, segment -> [taxon, segment, f] }
        .set{ch_input_meta}

    ch_input_fa
        .join(ch_input_meta, by: [0,1])
        .combine(ch_input.map{[it.taxon, it.exclusions]}, by: 0)
        .map{ taxon, segment, assembly, metadata, exclusions -> [taxon: taxon, segment: segment, assembly: assembly ? assembly : [], metadata: metadata ? metadata : [], exclusions: exclusions ? exclusions : []] }
        .set{ ch_input }

    ch_input.view()


    /*
    =============================================================================================================================
        CREATE REFERENCES
    =============================================================================================================================
    */
    if(params.create){

        CREATE_SUBWF(
            ch_input
        )
        ch_versions = ch_versions.mix(CREATE_SUBWF.out.versions)
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
