process mafft {
    label 'mafft'
    publishDir "${params.publish_path}/mafft", mode: params.publish_dir_mode

    input:
    path gene_family

    output:
    path "${gene_family.simpleName}.msa"

    script:
    """
    mafft --thread $task.cpus --auto --preservecase $gene_family > ${gene_family.simpleName}.msa
    """
}