process restore_duplicate_sequences {
    publishDir "${params.publish_path}/restore_duplicate_sequences", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa), path(dups_json)

    output:
    tuple val(index), path("${msa.name}.dups")

    script:
    """
    restore_duplicate_sequences.py -i $msa -j $dups_json -o ${msa.name}.dups
    """
}