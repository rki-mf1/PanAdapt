process codeml{
    publishDir "${params.publish_path}/codeml", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa), path(newick_tree), val(codeml_model)

    output:
    tuple val(index), path("${msa.simpleName}.txt"), val(codeml_model)
    tuple val(index), path("rst"), val(codeml_model)

    script:
    """
    sed 's|stewart\\.aa|${msa}|g;
    s|stewart\\.trees|${newick_tree}|g;
    s|out_M0\\.txt|${msa.baseName}.txt|g;
    s|NSsites_placeholder|${codeml_model}|g;' ${params.codeml_template} > codeml.ctl
    codeml codeml.ctl
    """
}