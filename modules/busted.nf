process busted {
    label 'hyphy'
    publishDir "${params.publish_path}/busted/", mode: params.publish_dir_mode
    
    input:
    tuple val(index), path(msa), path(newick_tree)

    output:
    tuple val(index), path("${index}.json")

    script:
    """
    hyphy busted --alignment $msa --tree $newick_tree --output ${index}.json
    """
}