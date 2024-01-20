process mask_ambiguities {
    label 'generic_small'
    publishDir "${params.publish_path}/mask_ambiguities", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa)

    output:
    tuple val(index), path("${msa.name}.masked")

    script:
    """
    mask_ambiguities.py -i $msa -o ${msa.name}.masked
    """
}