process RENAMEPAIREDFASTXHEADERS {
    tag "${meta.id}"
    label 'process_single'

    //conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.13'
        : 'quay.io/biocontainers/python:3.13'}"

    input:
    tuple val(meta), path(reads)
    path script

    output:
    tuple val(meta), path("*_renamed.*"), emit: reads
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def (in_f1, in_f2) = reads

    def prefix1 = in_f1.getName().tokenize('.')[0]
    def ext1 = in_f1.getName().tokenize('.')[1..-1].join('.')
    if (ext1[-3..-1]=='.gz') {
        ext1 = ext1.tokenize('.')[0..-2].join('.')
    }
    def out_f1 = "${prefix1}_renamed.${ext1}"

    def prefix2 = in_f2.getName().tokenize('.')[0]
    def ext2 = in_f2.getName().tokenize('.')[1..-1].join('.')
    if (ext2[-3..-1]=='.gz') {
        ext2 = ext2.tokenize('.')[0..-2].join('.')
    }
    def out_f2 = "${prefix2}_renamed.${ext2}"
    """
    python ${script} ${args} --no-zip-output -f ${in_f1} -r ${in_f2} -o ${out_f1} -l ${out_f2}
    gzip ${out_f1} ${out_f2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        renamepairedfastxheaders: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    def (in_f1, in_f2) = reads

    def prefix1 = in_f1.tokenize('.')[0]
    def ext1 = in_f1.tokenize('.')[1..-1].join('.')
    def out_f1 = "${prefix1}_renamed${ext1}"

    def prefix2 = in_f2.tokenize('.')[0]
    def ext2 = in_f2.tokenize('.')[1..-1].join('.')
    def out_f2 = "${prefix2}_renamed${ext2}"

    """
    touch ${out_f1}
    touch ${out_f2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        renamepairedfastxheaders: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
