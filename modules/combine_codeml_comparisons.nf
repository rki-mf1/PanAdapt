process combine_codeml_comparisons{
    publishDir "${params.publish_path}/combine_codeml_comparisons", mode: params.publish_dir_mode

    input:
    val model_comparisons

    output:
    path "model_comparison.tsv"

    script:
    def model_comparisons_str = model_comparisons.join(' ')
    """
    echo -e "Gene\tMO_vs_M1_p_value\tM1_vs_M2_p_value\tM7_vs_M8_p_value\tPreferred_Basic_Model\tPreferred_Site_Model" > model_comparison.tsv
    cat $model_comparisons_str >> model_comparison.tsv
    """
}