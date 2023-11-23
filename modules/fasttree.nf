process fasttree {
    publishDir "${params.publish_path}/fasttree", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa)

    output:
    tuple val(index), path("${msa.simpleName}.nwk")

    script:
    """
    fasttree -nt $msa > ${msa.simpleName}.nwk
    """
}