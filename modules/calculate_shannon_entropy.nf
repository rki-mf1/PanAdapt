process calculate_shannon_entropy {
    publishDir "${params.publish_path}/calculate_shannon_entropy/", mode: params.publish_dir_mode
    
    input:
    tuple val(index), path(site_table)

    output:
    tuple val(index), path("${index}.shannon.tsv")

    script:
    """
    calculate_shannon_entropy.py -i $site_table -o ${index}.shannon.tsv
    """
}