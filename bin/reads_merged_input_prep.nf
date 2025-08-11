import java.nio.file.Paths

def reads_merged_input_prep( reads_qc, cutadapt_channel ) {

    def dada2_input = reads_qc
        .mix(cutadapt_channel)                           // Combine fastp + cutadapt reads
        .map { meta, reads -> 
            [ groupKey(meta.subMap('id', 'single_end'), 2), meta, reads ]
        }
        .groupTuple(by: 0)                               // Group per sample (fastp + cutadapt)
        .map { key, metas, reads_list ->

            def meta = metas[0]                          // Take first meta

            //rename fastp_reads file names replacing the fastp extension
            // otherwise the READS_QC_MERGE_BEFOREHMM:FASTP fails because the input output name is the same
            def fastp_reads = reads_list[0] instanceof List
                ? reads_list[0].collect { file ->
                    def path = file instanceof Path ? file : Paths.get(file.toString())
                    def name = path.getFileName().toString().replaceFirst(/\.fastp/, '')
                    path.getParent().resolve(name)
                }
                : [ // for single-end reads, not a list
                    {
                        def path = reads_list[0] instanceof Path ? reads_list[0] : Paths.get(reads_list[0].toString())
                        def name = path.getFileName().toString().replaceFirst(/\.fastp/, '')
                        path.getParent().resolve(name)
                    }()
                ]
            
            def cutadapt_reads = reads_list[1] 

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
