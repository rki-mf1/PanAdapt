process msa_to_tsv {
    publishDir "${params.publish_path}/msa_to_tsv", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa)

    output:
    tuple val(index), path("${msa.simpleName}.tsv")

    script:
    """
    msa_to_tsv.py -i $msa -o ${msa.simpleName}.tsv
    """
}