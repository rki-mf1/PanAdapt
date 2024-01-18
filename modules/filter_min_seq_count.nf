process filter_min_seq_count {
    publishDir "${params.publish_path}/filter_min_seq_count", mode: params.publish_dir_mode

    input:
    path gene_family

    output:
    path "${gene_family.simpleName}.*"

    script:
    """
    seq_count=\$(grep -c '>' $gene_family || true)
    if [ "\$seq_count" -gt 2 ]; then
        cp $gene_family ${gene_family.simpleName}.filtered
    else
        cp $gene_family ${gene_family.simpleName}.tmp
    fi
    """
}