def reads_merged_input_prep( reads_qc, cutadapt_channel ) {

    def dada2_input = reads_qc
        .mix(cutadapt_channel)                           // Combine fastp + cutadapt reads
        .map { meta, reads -> 
            [ groupKey(meta.subMap('id', 'single_end'), 2), meta, reads ]
        }
        .groupTuple(by: 0)                               // Group per sample (fastp + cutadapt)
        .map { key, metas, reads_list ->

            def meta = metas[0]                          // Take first meta
            def fastp_reads = reads_list[0]              // FASTP output
            def cutadapt_reads = reads_list[1]           // Cutadapt output

            // Detect if cutadapt reads exist (SE vs PE check)
            def cutadapt_read_size = meta.single_end 
                                        ? cutadapt_reads.size() 
                                        : cutadapt_reads[0].size()

            // Select reads: cutadapt > fastp
            def final_reads = (cutadapt_read_size > 0) 
                                ? cutadapt_reads 
                                : fastp_reads

            [ meta, final_reads ]
        }

    return dada2_input
}
