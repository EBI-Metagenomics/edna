/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT EBI-METAGENOMICS MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READS_QC                     } from '../subworkflows/ebi-metagenomics/reads_qc/main'
include { READS_QC as READS_QC_MERGE   } from '../subworkflows/ebi-metagenomics/reads_qc/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_RAW         } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_CLEAN       } from '../modules/nf-core/fastqc/main'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../subworkflows/local/utils_nfcore_edna_pipeline'
include { PRIMER_IDENTIFICATION        } from '../subworkflows/local/primer_identification_swf.nf'
//include { AUTOMATIC_PRIMER_PREDICTION  } from '../subworkflows/local/automatic_primer_prediction.nf'
//include { CONCAT_PRIMER_CUTADAPT       } from '../subworkflows/local/concat_primer_cutadapt.nf'
include { MULTIQC                      } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EDNA {

    take:
    samplesheet // channel: samplesheet read in from --input
    main:
    
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Initialiase standard primer library for PIMENTO if user-given//
    // If there are no primers provided, it will fallback to use the default PIMENTO standard primer library
    std_primer_library = []

    if (params.std_primer_library){
        std_primer_library = file(params.std_primer_library, type: 'dir', checkIfExists: true)
    }

    FASTQC_RAW(
        samplesheet
    )
    fastqc_raw_html = FASTQC_RAW.out.html
    fastqc_raw_zip = FASTQC_RAW.out.zip
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())


    // Sanity checking and quality control of reads //
    READS_QC_MERGE(
        true, 
        samplesheet,
        true // merge
    )
    ch_versions = ch_versions.mix(READS_QC_MERGE.out.versions)

    // Run it again without merging to keep PE files unmerged for primer trimming+DADA2 //
    READS_QC(
        false, 
        samplesheet,
        false // merge
    )
    ch_versions = ch_versions.mix(READS_QC.out.versions)

    // Removes reads that passed sanity checks but are empty after QC with fastp //
    READS_QC_MERGE.out.reads_fasta.branch{ _meta, reads ->
                                qc_pass: reads.countFasta() > 0
                                qc_empty: reads.countFasta() == 0
                            }
                            .set { extended_reads_qc }
    
    FASTQC_CLEAN(
        READS_QC_MERGE.out.reads
    )
    fastqc_clean_html = FASTQC_CLEAN.out.html
    fastqc_clean_zip = FASTQC_CLEAN.out.zip
    ch_versions = ch_versions.mix(FASTQC_CLEAN.out.versions.first())

    // Identify whether primers exist or not in reads, separated by different amplified regions if more than one exists in a run //
    PRIMER_IDENTIFICATION(
        extended_reads_qc.qc_pass,
        std_primer_library
    )
    ch_versions = ch_versions.mix(PRIMER_IDENTIFICATION.out.versions)
     
    /*
    // Join primer identification flags with reads belonging to each run+amp_region //
    auto_trimming_input = PRIMER_IDENTIFICATION.out.conductor_out
                          .join(AMP_REGION_INFERENCE.out.extracted_var_out, by: [0])

    
    //Run subworkflow for automatic primer prediction
    //Outputs empty fasta file if no primers, or fasta file containing predicted primers
    
    AUTOMATIC_PRIMER_PREDICTION(
        auto_trimming_input
    )
    ch_versions = ch_versions.mix(AUTOMATIC_PRIMER_PREDICTION.out.versions)

    // Concatenate the different combinations of stranded std/auto primers for each run+amp_region //
    concat_input = PRIMER_IDENTIFICATION.out.std_primer_out
                   .join(AUTOMATIC_PRIMER_PREDICTION.out.auto_primer_trimming_out, by: [0])

    // Concatenate all primers for for a run, send them to cutadapt with original QCd reads for primer trimming //
    CONCAT_PRIMER_CUTADAPT(
        concat_input,
        READS_QC.out.reads
    )
    ch_versions = ch_versions.mix(CONCAT_PRIMER_CUTADAPT.out.versions)


    // Run the large dada2 imput preparation function //
    cutadapt_channel = CONCAT_PRIMER_CUTADAPT.out.cutadapt_out
                       .map { meta, reads -> 
                         [ meta.subMap('id', 'single_end'), meta['var_region'], meta['var_regions_size'], reads ]
                       }
    */
    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CLEAN.out.zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(READS_QC_MERGE.out.fastp_summary_json.map { it[1] })
    ch_versions = ch_versions.mix(FASTQC_CLEAN.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'edna_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
