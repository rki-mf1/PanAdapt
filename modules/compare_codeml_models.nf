process compare_codeml_models{
    publishDir "${params.publish_path}/compare_model_likelihood/${index}", mode: params.publish_dir_mode

    input:
    tuple val(index), val(lnL), val(np), val(codeml_model)

    output:
    path ".command.out"

    script:
    def lnL_str = lnL.join(' ')
    def np_str = np.join(' ')
    def codeml_model_str = codeml_model.join(' ')

    """
    compare_codeml_models.py -l $lnL_str -p $np_str -m $codeml_model_str -g ${index}
    """
}