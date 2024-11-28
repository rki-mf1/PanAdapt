process build_gene_families {
    label 'python'
    publishDir "${params.publish_path}/build_gene_families/", mode: params.publish_dir_mode

    input:
    path pan_matrix
    path fasta

    output:
    path "*.fna"
    path "*.faa"
    path "ambig_filter_counts.tsv"
    path "stop_filter_counts.tsv"

    script:
    """
    build_gene_families.py -m $pan_matrix -f $fasta -a 0.1
    """
}