
include { MAPSEQ                } from '../../modules/ebi-metagenomics/mapseq/main'
include { MAPSEQ2ASVTABLE       } from '../../modules/local/mapseq2asvtable/main.nf'
include { MAKE_ASV_COUNT_TABLES } from '../../modules/local/make_asv_count_tables/main.nf'
include { KRONA_KTIMPORTTEXT    } from '../../modules/ebi-metagenomics/krona/ktimporttext/main'
include { EXTRACT_ASVS_LEFT     } from '../../modules/local/extract_asvs_left/main.nf'

workflow MAPSEQ_ASV_KRONA {
    
    take:
        dada2_output
        krona_tuple

    main:

        ch_versions = Channel.empty()

        mapseq_input = dada2_output
                       .map { meta, maps, asv_seqs, filt_reads ->
                            [ meta, asv_seqs ]
                        }
        MAPSEQ(
            mapseq_input,
            tuple(krona_tuple[0], krona_tuple[1], krona_tuple[3])
        )
        ch_versions = ch_versions.mix(MAPSEQ.out.versions.first())


        MAPSEQ2ASVTABLE(
            MAPSEQ.out.mseq,
            krona_tuple[4] // db_label
        )
        ch_versions = ch_versions.mix(MAPSEQ2ASVTABLE.out.versions.first())

        /*
        // Transpose by var region in case any samples have more than one
        split_mapseq2asvtable = MAPSEQ2ASVTABLE.out.asvtaxtable
                                .map { meta, asvtaxtable -> 
                                     [ meta.subMap('id', 'single_end'), asvtaxtable ]       
                                    }

        // Transpose by var region in case any samples have more than one. Also reorder the inputs slightly
        split_input = dada2_output
            .map { meta, maps, asv_seqs, filt_reads -> 
                [ meta.subMap('id', 'single_end'), maps, filt_reads ]
            }
            .join(split_mapseq2asvtable, by: 0)
            .map { meta, maps, filt_reads, asvtaxtable ->
                [ meta, maps, asvtaxtable, filt_reads ] 
            }

        */
        // Add in the concatenated var region channel to the rest of the input
        final_asv_count_table_input = dada2_output
            .join(MAPSEQ2ASVTABLE.out.asvtaxtable, by: 0)
            .map { meta, maps, asv_seqs, filt_reads, asvtaxtable ->
                tuple( meta + [ 'var_region': 'all' ], maps, asvtaxtable, filt_reads, asv_seqs )
            }

        
        MAKE_ASV_COUNT_TABLES(
            final_asv_count_table_input,
            krona_tuple[4]
        )
        ch_versions = ch_versions.mix(MAKE_ASV_COUNT_TABLES.out.versions.first())

/*
        KRONA_KTIMPORTTEXT(
            MAKE_ASV_COUNT_TABLES.out.asv_count_tables_out,
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions.first())
   */
    emit:
      //  asv_count_tables_out = MAKE_ASV_COUNT_TABLES.out.asv_count_tables_out
       // asvs_left = MAKE_ASV_COUNT_TABLES.out.asv_read_counts_out
     //   asvtaxtable = MAPSEQ2ASVTABLE.out.asvtaxtable
     //   krona_out = KRONA_KTIMPORTTEXT.out.html
        versions = ch_versions
}
