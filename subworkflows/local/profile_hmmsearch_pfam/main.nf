include { FASTAEMBEDLENGTH } from '../../../modules/local/fastaembedlength/main'
include { SEQKIT_TRANSLATE } from '../../../modules/nf-core/translate/main'
include { HMMER_HMMSEARCH } from '../../../modules/nf-core/hmmsearch/main'
include { PARSEHMMSEARCHCOVERAGE } from '../../../modules/local/parsehmmsearchcoverage/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'

workflow PROFILE_HMMSEARCH_PFAM {

    take:
    reads_fasta
    pfam_db
    reads_json

    main:
    ch_versions = Channel.empty()
    FASTAEMBEDLENGTH(reads_fasta, file("${projectDir}/bin/fastx_embed_length.py"))
    
    SEQKIT_TRANSLATE(FASTAEMBEDLENGTH.out.fasta)
    
    /*
    * Adaptive FASTA chunking for HMMER parallelization
    * 
    * Uses ideal chunk size (better performance) for small files,
    * but switches to larger chunks for big files to respect the
    * maximum chunk limit (prevents job scheduler overload).
    * 
    * Sequence-based splitting ensures predictable processing times
    * and prevents truncated sequences.
    */

    ch_chunked_pfam_in = SEQKIT_TRANSLATE.out.fastx
        .flatMap{ meta, fasta ->
            def totalSeqs = fasta.countFasta()
            def maxChunks = params.hmmsearch_max_chunks ?: 200
            def idealSeqsPerChunk = params.hmmsearch_seqs_per_chunk ?: 1000

            // Adaptive strategy: ideal chunks for small files, larger chunks for big files
            def seqsPerChunk = (totalSeqs <= maxChunks * idealSeqsPerChunk) ? 
                idealSeqsPerChunk :                           
                Math.ceil(totalSeqs / maxChunks) as Integer   
                
            def chunks = fasta.splitFasta(by: seqsPerChunk, file: true)
            
            chunks.collect{ chunk -> tuple(groupKey(meta, chunks.size()), chunk) }
        }
        .combine(pfam_db)
        .map{ meta, reads, db -> [meta, db, reads, true, true, true] }

    HMMER_HMMSEARCH(ch_chunked_pfam_in)
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

    COMBINEHMMSEARCHTBL(
        HMMER_HMMSEARCH.out.domain_summary.groupTuple()
    )

    PARSEHMMSEARCHCOVERAGE(COMBINEHMMSEARCHTBL.out.concatenated_result.join(reads_json), file("${projectDir}/bin/hmmer_domtbl_parse_coverage.py"))
    ch_versions = ch_versions.mix(PARSEHMMSEARCHCOVERAGE.out.versions)

    emit:
    profile  = PARSEHMMSEARCHCOVERAGE.out.tsv
    domtbl   = COMBINEHMMSEARCHTBL.out.concatenated_result
    versions = ch_versions                                   // channel: [ versions.yml ]
}

