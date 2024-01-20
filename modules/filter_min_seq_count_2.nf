process filter_min_seq_count_2 {
    label 'generic_small'
    publishDir "${params.publish_path}/filter_min_seq_count_2", mode: params.publish_dir_mode

    input:
    tuple val(index), path(gene_family)

    output:
    tuple val(index), path("${gene_family.simpleName}.*")

    script:
    """
    seq_count=\$(grep -c '>' $gene_family)
    if [ "\$seq_count" -gt 2 ]; then
        cp $gene_family ${gene_family.simpleName}.filtered
    else
        cp $gene_family ${gene_family.simpleName}.tmp
    fi
    """
}