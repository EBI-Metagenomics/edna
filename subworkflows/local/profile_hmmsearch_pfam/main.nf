include { FASTAEMBEDLENGTH } from '../../../modules/local/fastaembedlength/main'
include { SEQKIT_TRANSLATE } from '../../../modules/nf-core/translate/main'
include { HMMER_HMMSEARCH } from '../../../modules/nf-core/hmmsearch/main'
include { PARSEHMMSEARCHCOVERAGE } from '../../../modules/local/parsehmmsearchcoverage/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'

workflow PROFILE_HMMSEARCH_PFAM {

    take:
    reads_fasta
    pfam_db

    main:
    ch_versions = Channel.empty()
    FASTAEMBEDLENGTH(reads_fasta, file("${projectDir}/bin/fastx_embed_length.py"))
    SEQKIT_TRANSLATE(FASTAEMBEDLENGTH.out.fasta)

    ch_chunked_pfam_in = SEQKIT_TRANSLATE.out.fastx
        .flatMap{ meta, fasta ->
             def chunks = fasta.splitFasta(file: true, size: params.hmmsearch_chunksize)
             chunks.collect{ chunk -> tuple(groupKey(meta, chunks.size()), chunk) }
         }
        .combine(pfam_db)
        .map{ meta, reads, db -> [meta, db, reads, true, true, true] }

    HMMER_HMMSEARCH(ch_chunked_pfam_in)
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

    COMBINEHMMSEARCHTBL(
        HMMER_HMMSEARCH.out.domain_summary.groupTuple()
    )

    PARSEHMMSEARCHCOVERAGE(COMBINEHMMSEARCHTBL.out.concatenated_result, file("${projectDir}/bin/hmmer_domtbl_parse_coverage.py"))
    ch_versions = ch_versions.mix(PARSEHMMSEARCHCOVERAGE.out.versions)

    emit:
    profile  = PARSEHMMSEARCHCOVERAGE.out.tsv
    versions = ch_versions                     // channel: [ versions.yml ]
}

