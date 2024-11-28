process merge_fubar_stats {
    publishDir "${params.publish_path}/merge_fubar_stats", mode: params.publish_dir_mode

    input:
    path fubar_jsons

    output:
    file "merged_fubar_stats.tsv"

    script:
    """
    printf '%s\\n' ${fubar_jsons} > file_list.tmp
    first_file=\$(head -n 1 file_list.tmp)
    head -n 1 "\$first_file" > merged_fubar_stats.tsv
    while read file; do
        tail -n 1 "\$file" >> merged_fubar_stats.tsv
    done < file_list.tmp
    """
}