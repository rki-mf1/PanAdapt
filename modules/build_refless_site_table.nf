process build_refless_site_table {
    label 'python'
    publishDir "${params.publish_path}/build_refless_site_table", mode: params.publish_dir_mode
    scratch true
    
    input:
    tuple val(index), path(msa)

    output:
    tuple val(index), path("${index}.tsv")

    script:
    """
    build_site_table.py -i $msa -o ${index}.tsv
    """
}