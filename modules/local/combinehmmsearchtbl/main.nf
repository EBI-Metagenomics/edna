process COMBINEHMMSEARCHTBL {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(input_files, stageAs: "inputs/?/*")

    output:
    tuple val(meta), path("${output_name}"), emit: concatenated_result

    script:
    def args = task.ext.args ?: ''
    output_name = (input_files instanceof List ? input_files[0] : input_files).name.tokenize('/')[-1]

    """
    export inputs=\$(find -L inputs -name "*tbl.gz")
    export first_file=\$(echo \$inputs | head -n 1)
    gunzip -c \${first_file} | grep '#' > header.txt || true
    cat \$inputs | gunzip | grep -v '#' > without_header.txt || true
    cat header.txt without_header.txt | gzip > ${output_name} || true
    """

    stub:
    def args = task.ext.args ?: ''
    output_name = input_files[0].name.tokenize('/')[-1]

    """
    touch ${output_name}
    """
}
