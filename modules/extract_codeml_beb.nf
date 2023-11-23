process extract_codeml_beb {
    publishDir "${params.publish_path}/extract_codeml_beb", mode: params.publish_dir_mode

    input:
    tuple val(index), path(rst), val(codeml_model)

    output:
    tuple val(index), path("${index}.*")

    script:
    """
    if [ "$codeml_model" == '8' ]; then 
        extract_codeml_beb.py -i $rst -o ${index}.beb
    else 
        touch ${index}.tmp
    fi
    """
}