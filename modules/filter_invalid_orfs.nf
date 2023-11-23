process filter_invalid_orfs {
    publishDir "${params.publish_path}/filter_invalid_orfs/${gff.baseName}", mode: params.publish_dir_mode

    input:
    path gff

    output:
    path "${gff.baseName}_filtered.gff"

    script:
    """
    filter_invalid_orfs.py -i $gff -o ${gff.baseName}_filtered.gff
    """
}