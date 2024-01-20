process build_site_fubar_table {
    label 'generic_small'
    publishDir "${params.publish_path}/build_site_fubar_table/", mode: params.publish_dir_mode
    
    input:
    tuple val(index), path(site_table), path(fubar_json)

    output:
    tuple val(index), path("${index}.fubar.tsv")

    script:
    """
    build_site_fubar_table.py -i $site_table -j $fubar_json -o ${index}.fubar.tsv
    """
}