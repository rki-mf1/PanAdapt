process reformat_matrix {
    label 'python'
    publishDir "${params.publish_path}/reformat_matrix/", mode: params.publish_dir_mode

    input:
    path pan_matrix

    output:
    path "reformatted_matrix.tsv"

    script:
    """
    reformat_matrix.py -i $pan_matrix -o "reformatted_matrix.tsv"
    """
    }