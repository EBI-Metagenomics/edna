process FASTAEMBEDLENGTH {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.13'
        : 'quay.io/biocontainers/python:3.13'}"

    input:
    tuple val(meta), path(fasta)
    path script

    output:
    tuple val(meta), path("${out_fn}"), emit: fasta
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def base = fasta.getName().tokenize('.')[0..-2].join('.')
    def ext = fasta.getName().tokenize('.')[-1]
    out_fn = "${base}.renamed.${ext}"

    """
    python ${script} -i ${fasta} -o ${out_fn} --output_gzip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def base = fasta.getName().tokenize('.')[0..-2].join('.')
    def ext = fasta.getName().tokenize('.')[-1]
    out_fn = "${base}.renamed.${ext}"
    
    """
    touch ${out_fn}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
