process extract_fubar {
    label 'python'
    publishDir "${params.publish_path}/extract_fubar/", mode: params.publish_dir_mode

    input:
    tuple val(index), path(fubar_json)

    output:
    tuple val(index), path("${index}.tsv")
    path "${index}_fubar_stats.tsv"

    script:
    """
    extract_fubar.py -i $fubar_json -g $index -o ${index}.tsv -s ${index}_fubar_stats.tsv
    """
}