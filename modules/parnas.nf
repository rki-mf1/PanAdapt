process parnas {
    publishDir "${params.output}/parnas", mode: params.publish_dir_mode

    input:
    tuple val(index), path(newick_tree)

    output:
    tuple val(index), path("${newick_tree.simpleName}.*")

    script:
    """
    seq_count=`grep -o ',' $newick_tree | wc -l`
    if [ "\$seq_count" -le "${params.parnas_threshold}" ]; then
        cp $newick_tree ${newick_tree.simpleName}.parnas
    else
        parnas -t $newick_tree --cover --radius 0.001 --subtree ${newick_tree.simpleName}.parnas_reduced
        sed -i "s/'//g" ${newick_tree.simpleName}.parnas_reduced
    fi
    """
}