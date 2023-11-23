process extract_codeml_parameters {
    publishDir "${params.publish_path}/extract_codeml_parameters", mode: params.publish_dir_mode

    input:
    tuple val(index), path(codeml_output), val(codeml_model)

    output:
    stdout

    script:
    """
    extract_codeml_parameters.py $codeml_output $index $codeml_model
    """
}
