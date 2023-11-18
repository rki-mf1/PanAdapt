process nexus_to_newick {
    publishDir "${params.output}/nexus_to_newick", mode: params.publish_dir_mode

    input:
    tuple val(index), path(newick_tree)

    output:
    tuple val (index), path("${newick_tree.simpleName}.nwk_reduced")

    script:
    """
    if head -n 1 $newick_tree | grep -q '^#NEXUS'; then
        nexus_to_newick.py -i $newick_tree -o ${newick_tree.simpleName}.nwk_reduced
    else
        cp $newick_tree ${newick_tree.simpleName}.nwk_reduced
    fi
    """
}