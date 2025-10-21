
include { STD_PRIMER_FLAG     } from '../../modules/local/std_primer_flag/main.nf'


workflow PRIMER_IDENTIFICATION {
    
    // Subworkflow that attempts to identify the presence of primers in sequencing files (fastq)
    // Checks for the presence of standard primers from a known library using params.std_primer_library

    take:
        reads_merged
        std_primer_library
    main:

        ch_versions = Channel.empty()

        // Standard Library primers
        STD_PRIMER_FLAG(
            reads_merged,
            std_primer_library
        )
        ch_versions = ch_versions.mix(STD_PRIMER_FLAG.out.versions.first())
        

    emit:
        std_primer_out = STD_PRIMER_FLAG.out.std_primer_out
        versions = ch_versions
    
}