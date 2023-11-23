process remove_duplicate_sequences {
    publishDir "${params.publish_path}/remove_duplicate_sequences", mode: params.publish_dir_mode

    input:
    path gene_family

    output:
    path "${gene_family.name}.no_dups"
    path "${gene_family.name}.json"

    script:
    """
    remove_duplicate_sequences.py -i $gene_family -o ${gene_family.name}.no_dups -j ${gene_family.name}.json -r ${params.ref_id} 
    """
}
