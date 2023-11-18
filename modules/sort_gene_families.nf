process sort_gene_families {
    publishDir "${params.output}/sort_gene_families", mode: params.publish_dir_mode
    input:
    path gene_family

    output:
    path "${gene_family.name}.sorted" 

    script:
    """
    seqkit sort $gene_family > ${gene_family.name}.sorted
    """
}