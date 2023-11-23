process reduce_msa {
    publishDir "${params.publish_path}/reduce_msa", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa), path(newick_tree)

    output:
    tuple val(index), path("${msa.simpleName}.reduced_msa")

    script:
    """
    reduce_msa.py -i $msa -t $newick_tree -o ${msa.simpleName}.reduced_msa
    """
}