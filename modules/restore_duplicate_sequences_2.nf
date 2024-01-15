process restore_duplicate_sequences_2 {
    publishDir "${params.publish_path}/restore_duplicate_sequences_2", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa), path(dups_json)

    output:
    tuple val(index), path("${msa.name}.dups_2")

    script:
    """
    restore_duplicate_sequences.py -i $msa -j $dups_json -o ${msa.name}.dups_2
    """
}