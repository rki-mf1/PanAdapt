process build_site_stats {
    label 'python'
    publishDir "${params.publish_path}/build_site_stats/", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa), path(fubar)

    output:
    tuple val(index), path("${index}_site_stats.tsv")

    script:
    """
    build_site_stats.py -i $msa -f $fubar -o ${index}_site_stats.tsv
    """
}