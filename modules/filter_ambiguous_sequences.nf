process filter_ambiguous_sequences {
    publishDir "${params.publish_path}/filter_ambiguous_sequences", mode: params.publish_dir_mode

    input:
    path gene_family

    output:
    path "${gene_family.name}.no_ambig"

    script:
    """
    filter_ambiguous_sequences.py -i $gene_family -o ${gene_family.name}.no_ambig
    """
}