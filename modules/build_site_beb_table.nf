process build_site_beb_table {
    publishDir "${params.publish_path}/build_site_beb_table", mode: params.publish_dir_mode

    input:
    tuple val(index), path(per_site_table), path(beb_table)

    output:
    tuple val(index), path("${per_site_table.name}.beb")

    script:
    """
    build_site_beb_table.py -i $per_site_table -b $beb_table -o ${per_site_table.name}.beb
    """
}