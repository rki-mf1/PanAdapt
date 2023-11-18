process translate_gene_families {
    publishDir "${params.output}/translate_gene_families", mode: params.publish_dir_mode

    input:
    path gene_family

    output:
    path "${gene_family.name}.translated"

    script:
    """
    transeq -sequence $gene_family -outseq ${gene_family.name}.translated
    sed -i 's/_1\$//' ${gene_family.name}.translated
    """
}