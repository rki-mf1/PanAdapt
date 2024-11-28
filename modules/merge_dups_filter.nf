process merge_dups_filter {
    publishDir "${params.publish_path}/merge_dups_filter", mode: params.publish_dir_mode

    input:
    path dups_filters

    output:
    file "merged_dups_filters.tsv"

    script:
    """
    printf '%s\\n' ${dups_filters} > file_list.tmp
    first_file=\$(head -n 1 file_list.tmp)
    head -n 1 "\$first_file" > merged_dups_filters.tsv
    while read file; do
        tail -n 1 "\$file" >> merged_dups_filters.tsv
    done < file_list.tmp
    """
}