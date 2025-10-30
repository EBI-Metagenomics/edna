/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT EBI-METAGENOMICS MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READS_QC                     } from '../subworkflows/ebi-metagenomics/reads_qc/main'
include { READS_QC as READS_QC_MERGE   } from '../subworkflows/ebi-metagenomics/reads_qc/main'
include { READS_QC as READS_QC_MERGE_BEFOREHMM   } from '../subworkflows/ebi-metagenomics/reads_qc/main'


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
include { CONCAT_PRIMER_CUTADAPT       } from '../subworkflows/local/concat_primer_cutadapt.nf'
include { PROFILE_HMMSEARCH_PFAM       } from '../subworkflows/local/profile_hmmsearch_pfam/main'
include { DADA2_SWF                    } from '../subworkflows/local/dada2_swf.nf'
include { MAPSEQ_ASV_KRONA             } from '../subworkflows/local/mapseq_asv_krona_swf.nf'
include { MULTIQC                      } from '../modules/nf-core/multiqc/main'

// Import samplesheetToList from nf-schema //
include { samplesheetToList            } from 'plugin/nf-schema'

// Import reads_merged_input_prep function (it's very big and deserved to be in its own file) //
include { reads_merged_input_prep      } from '../bin/reads_merged_input_prep.nf'

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

     
     /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INITIALISE REFERENCE DATABASE INPUT TUPLES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Regular ASV resolution method //
    dada2_krona_bold_tuple = tuple(
        file(params.bold_db_fasta, checkIfExists: true),
        file(params.bold_db_tax, checkIfExists: true),
        file(params.bold_db_otu, checkIfExists: true),
        file(params.bold_db_mscluster, checkIfExists: true),
        params.dada2_bold_label
    )

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
        false, 
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

    // Filter and branch reads based on minimum read count with logging
    READS_QC_MERGE.out.reads_se_and_merged.branch{ meta, reads ->
                                    def count = reads.countFastq()
                                    qc_pass: count >= params.min_read_count
                                    qc_fail: count < params.min_read_count
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
     
    // Concatenate all primers for for a run, send them to cutadapt with original QCd reads for primer trimming //
    CONCAT_PRIMER_CUTADAPT(
        PRIMER_IDENTIFICATION.out.std_primer_out,
        READS_QC.out.reads
    )
    ch_versions = ch_versions.mix(CONCAT_PRIMER_CUTADAPT.out.versions)
     
    reads_merge_input = reads_merged_input_prep(READS_QC.out.reads, CONCAT_PRIMER_CUTADAPT.out.cutadapt_out)
    
    READS_QC_MERGE_BEFOREHMM(
        false, 
        reads_merge_input,
        true // merge
    )
    ch_versions = ch_versions.mix(READS_QC_MERGE_BEFOREHMM.out.versions) 

    // Pfam profiling
    pfam_db = params.pfam_coi_db ?
    Channel
        .fromPath(params.pfam_coi_db, checkIfExists: true)
        .first() :
    Channel.empty()
    
    PROFILE_HMMSEARCH_PFAM(
        READS_QC_MERGE_BEFOREHMM.out.reads_fasta,
        pfam_db,
        READS_QC_MERGE_BEFOREHMM.out.fastp_summary_json
    )
    ch_versions = ch_versions.mix(PROFILE_HMMSEARCH_PFAM.out.versions)
        
    // Filter samples based on reads_percentage threshold and get filtered domtbl
    ch_passed_samples = PROFILE_HMMSEARCH_PFAM.out.profile
        .filter { meta, tsv_file ->
            def threshold = params.reads_percentage_threshold ?: 0.70
            
            try {
                def lines = tsv_file.readLines()
                def dataLine = lines[1] // Skip header, get the single data row
                def columns = dataLine.split('\t')
                def readsPercentageStr = columns[4]
                
                def readsPercentage = readsPercentageStr as Double
                def passes = readsPercentage >= threshold
                
                return passes
                
            } catch (Exception e) {
                return false
            }
        }
        .map { meta, tsv_file -> meta }

    ch_filtered_domtbl = ch_passed_samples.join(PROFILE_HMMSEARCH_PFAM.out.domtbl) 

    // Run DADA2 ASV generation with filtered samples //
    DADA2_SWF(
        reads_merge_input,
        ch_filtered_domtbl
    )
    ch_versions = ch_versions.mix(DADA2_SWF.out.versions)
 
    def dada2_stats_fail = DADA2_SWF.out.dada2_stats_fail.map { meta, stats_fail ->
                                def key = meta.subMap('id', 'single_end')
                                return [key, stats_fail]
                            }

    // ASV taxonomic assignments + generate Krona plots for each run+amp_region //
    MAPSEQ_ASV_KRONA(
        DADA2_SWF.out.dada2_out,
        dada2_krona_bold_tuple
    )
    ch_versions = ch_versions.mix(MAPSEQ_ASV_KRONA.out.versions)

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

    /*****************************/
    /* End of execution reports */
    /****************************/
 
    // Extract runs that failed SeqFu check //
    READS_QC.out.seqfu_check
        .splitCsv(sep: "\t", elem: 1)
        .filter { meta, seqfu_res ->
            seqfu_res[0] != "OK"
        }
        .map { meta, __ -> "${meta.id},seqfu_fail" }
        .set { seqfu_fails }

    // Extract runs that failed Suffix Header check //
    READS_QC.out.suffix_header_check
        .filter { meta, sfxhd_res ->
            sfxhd_res.countLines() != 0
        }
        .map { meta, __ -> "${meta.id},sfxhd_fail"  }
        .set { sfxhd_fails }

    // Extract runs that failed Library Strategy check //
    READS_QC_MERGE.out.amplicon_check
        .filter { meta, strategy ->
            strategy != "AMPLICON"
        }
        .map { meta, __ -> "${meta.id},libstrat_fail" }
        .set { libstrat_fails }

    // Extract runs that had zero reads after fastp //
    extended_reads_qc.qc_fail.map { meta, __ -> "${meta.id},min_reads"  }
        .set { min_reads_fails }

    // Extract runs that failed the reads_percentage_threshold parameter
    PROFILE_HMMSEARCH_PFAM.out.profile
        .map { meta, tsv_file -> [meta.id, meta] }
        .join(ch_passed_samples.map { meta -> [meta.id, meta] }, remainder: true)
        .filter { id, meta, passed_meta ->
            passed_meta == null  // These are samples that didn't pass
        }
        .map { id, meta, __ -> "${meta.id},reads_percentage_fail" }
        .set { reads_percentage_fails }

    // Save all failed runs to file //
    all_failed_runs = seqfu_fails.concat( sfxhd_fails, libstrat_fails, min_reads_fails, reads_percentage_fails)
    all_failed_runs.collectFile(name: "qc_failed_runs.csv", storeDir: "${params.outdir}", newLine: true, cache: false)

    // Extract passed runs, describe whether those passed runs also ASV results //
    DADA2_SWF.out.dada2_report.map { meta, dada2_report -> [ ["id": meta.id, "single_end": meta.single_end], dada2_report ] }
    .concat(
        ch_passed_samples.map { meta -> [["id": meta.id, "single_end": meta.single_end], "qc_pass"] },
        dada2_stats_fail
    )
    .groupTuple()
    .map { meta, results ->
        if ( results.size() == 3 ) {
            return "${meta.id},all_results"
        }
        else {
            if (results.find { it == "true" }) {
                return "${meta.id},dada2_stats_fail"
            } else {
                return "${meta.id},no_asvs"
            }
        }
        error "Unexpected. meta: ${meta}, results: ${results}"
    }
    .set { final_passed_runs }

    // Save all passed runs to file //
    final_passed_runs.collectFile(name: "qc_passed_runs.csv", storeDir: "${params.outdir}", newLine: true, cache: false)
    .set { passed_runs_path }


    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
