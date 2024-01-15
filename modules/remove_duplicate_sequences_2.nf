process remove_duplicate_sequences_2 {
    publishDir "${params.publish_path}/remove_duplicate_sequences_2", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa)

    output:
    tuple val(index), path("${msa.name}.no_dups")
    tuple val(index), path("${msa.name}.json")

    script:
    """
    remove_duplicate_sequences.py -i $msa -o ${msa.name}.no_dups -j ${msa.name}.json -r ${params.ref_id} 
    """
}