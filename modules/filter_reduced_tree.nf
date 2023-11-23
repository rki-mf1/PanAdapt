process filter_reduced_tree{
    publishDir "${params.publish_path}/filter_reduced_tree", mode: params.publish_dir_mode

    input:
    tuple val(index), path(reduced_tree)

    output:
    tuple val(index), path("${reduced_tree.simpleName}.*")

    script:
    """
    set -xe
    seq_count=\$(grep -o ',' $reduced_tree | wc -l)
    if [ "\$seq_count" -gt 0 ]; then
        if grep -q "^#NEXUS" $reduced_tree; then
            grep "^    TREE 1 =" $reduced_tree | sed -n 's/^.*TREE 1 = //p' > ${reduced_tree.simpleName}.filtered
        else
            cp $reduced_tree ${reduced_tree.simpleName}.filtered
        fi
    else
        cp $reduced_tree ${reduced_tree.simpleName}.tmp
    fi
    """
}
