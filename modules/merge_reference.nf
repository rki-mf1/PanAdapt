process merge_reference {
    label 'python'
    publishDir "${params.publish_path}/merge_reference/", mode: params.publish_dir_mode

    input:
    tuple val(index), path(site_stats), path(aligned_reference)

    output: 
    tuple val(index), path("${index}_site_stats_ref.tsv")

    script:
    """
    merge_reference.py -i $site_stats -r $aligned_reference -o ${index}_site_stats_ref.tsv
    """
}