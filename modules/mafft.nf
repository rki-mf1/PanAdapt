process mafft {
    publishDir "${params.publish_path}/mafft", mode: params.publish_dir_mode

    input:
    path gene_family

    output:
    path "${gene_family.simpleName}.msa"

    script:
    """
    mafft --thread -1 --auto --preservecase $gene_family > ${gene_family.simpleName}.msa
    """
}