process merge_gene_stats {
    label 'python'
    publishDir "${params.publish_path}/merge_gene_stats/", mode: params.publish_dir_mode

    input:
    path reformatted_matrix
    path ambig_filter_counts
    path stop_filter_counts
    path merged_dups_filter
    path merged_fubar_stats

    output:
    path "merged_gene_stats.tsv"

    script:
    """
    merge_gene_stats.py -i $reformatted_matrix \
                        -a $ambig_filter_counts \
                        -s $stop_filter_counts \
                        -d $merged_dups_filter \
                        -f $merged_fubar_stats \
                        -o merged_gene_stats.tsv
    """
}