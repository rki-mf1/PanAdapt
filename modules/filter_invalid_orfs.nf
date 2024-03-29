process filter_invalid_orfs {
    label 'python'
    publishDir "${params.publish_path}/filter_invalid_orfs/${gff.baseName}", mode: params.publish_dir_mode

    input:
    path gff

    output:
    path "${gff.name}_filtered"

    script:
    """
    filter_invalid_orfs.py -i $gff -o ${gff.name}_filtered
    """
}