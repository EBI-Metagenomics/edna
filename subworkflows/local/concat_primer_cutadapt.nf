
include { REV_COMP_SE_PRIMERS     } from '../../modules/local/rev_comp_se_primers/main.nf'
include { SPLIT_PRIMERS_BY_STRAND } from '../../modules/local/split_primers_by_strand.nf'
include { CUTADAPT                } from '../../modules/ebi-metagenomics/cutadapt/main.nf'

workflow CONCAT_PRIMER_CUTADAPT {
    
    take:
        concat_input
        reads
    main:

        ch_versions = Channel.empty()
        
        REV_COMP_SE_PRIMERS(
            concat_input
        )
        ch_versions = ch_versions.mix(REV_COMP_SE_PRIMERS.out.versions.first())
    
        SPLIT_PRIMERS_BY_STRAND(
            REV_COMP_SE_PRIMERS.out.rev_comp_se_primers_out
        )
        ch_versions = ch_versions.mix(SPLIT_PRIMERS_BY_STRAND.out.versions.first())

        // Join concatenated primers to the fastp-cleaned paired reads files and run cutadapt on them
        cutadapt_input = SPLIT_PRIMERS_BY_STRAND.out.stranded_primer_out
                        .map{ meta, fwd_primer, rev_primer ->
                            [ meta.subMap('id', 'single_end'), meta['var_region'], meta['var_regions_size'], [fwd_primer, rev_primer] ]
                        }
                        .join(reads, by: [0])
                        .map{ meta, var_region, var_regions_size, primers, final_reads -> 
                            [ meta + ["var_region":var_region, "var_regions_size": var_regions_size], final_reads, primers ]
                        }

        CUTADAPT(
            cutadapt_input
        )
        ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

    emit:
        cutadapt_out = CUTADAPT.out.reads
        cutadapt_json = CUTADAPT.out.json
        versions = ch_versions

}