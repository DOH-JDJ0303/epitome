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
                segment:  it.containsKey('segment') ? ( it.segment ? it.segment : 'wg') : 'wg',
                assembly: it.containsKey('assembly') ? ( it.assembly ? file(it.assembly, checkIfExists: true) : [] ) : [], 
                metadata: it.containsKey('metadata') ? ( it.metadata ? file(it.metadata, checkIfExists: true) : [] ) : [],
                segmentSynomyms: it.containsKey('segmentSynonyms') ? ( it.segmentSynonyms ? it.segmentSynonyms : null ) : null,
                segmented: it.containsKey('segmented') ? ( it.segmented ? it.segmented : 'FALSE' ) : 'FALSE',
                exclusions: it.containsKey('exclusions') ? ( it.exclusions ? file(it.exclusions, checkIfExists: true) : [] ) : [],
                ]
                }
        .set{ ch_input_full }

        ch_input_full.map{ [ taxon: it.taxon, segment: it.segment, assembly: it.assembly, metadata: it.metadata ] }.set{ ch_input_part }

    /*
    =============================================================================================================================
        GATHER NCBI DATA
    =============================================================================================================================
    */
    if(params.ncbi){

        NCBI_DATA_SUBWF(
            ch_input_full
        )
        ch_versions = ch_versions.mix(NCBI_DATA_SUBWF.out.versions)
        
        NCBI_DATA_SUBWF
            .out
            .input
            .concat(ch_input_part.filter{ it.segment && it.assembly })
            .map{ [ it.taxon, it.segment, it.assembly, it.metadata ] }
            .groupTuple(by: [0,1])
            .set{ch_input_part}
    }

    /*
    =============================================================================================================================
        MERGE INPUTS
    =============================================================================================================================
    */
    // Reformat and merge inputs (from NCBI and the user)
    MERGE_INPUTS (
        ch_input_part
    )
    ch_versions = ch_versions.mix(MERGE_INPUTS.out.versions.first())
    MERGE_INPUTS
        .out
        .merged
        .join( ch_input_full.map{ [ it.taxon, it.segment, it.exclusions ] }, by: [0,1], remainder: true)
        .filter{ it[2] }
        .map{ taxon, segment, assembly, metadata, exclusions -> [ taxon: taxon, segment: segment, assembly: assembly, metadata: metadata, exclusions: exclusions ? exclusions : [] ] }
        .set{ ch_input }

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
