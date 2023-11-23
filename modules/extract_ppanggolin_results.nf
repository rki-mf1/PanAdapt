process extract_ppanggolin_results {
    publishDir "${params.publish_path}/extract_ppanggolin_results/", mode: params.publish_dir_mode

    input:
    path matrix_csv

    output:
    path "ppanggolin_results.tsv"

    script:
    """
    extract_ppanggolin_results.py -i $matrix_csv -o ppanggolin_results.tsv
    """
}
