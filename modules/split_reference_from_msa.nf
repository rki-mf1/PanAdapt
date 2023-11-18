process split_reference_from_msa{
    publishDir "${params.output}/split_reference_from_msa", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa)

    output:
    tuple val(index), path("${msa.name}.no_ref")
    tuple val(index), path("${msa.baseName}.ref")

    """
    split_reference_from_msa.py -i $msa -r ${params.ref_id} -o ${msa.name}.no_ref -u ${msa.baseName}.ref
    """
}