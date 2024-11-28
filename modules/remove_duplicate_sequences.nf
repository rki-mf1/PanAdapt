process remove_duplicate_sequences {
    label 'python'
    publishDir "${params.publish_path}/remove_duplicate_sequences", mode: params.publish_dir_mode

    input:
    tuple val(index), path(msa)

    output:
    tuple val(index), path("${msa.name}.no_dups") optional true
    path "${index}_dups_filter_counts.tsv"

    script:
    """
    remove_duplicate_sequences.py -i $msa -o ${msa.name}.no_dups -g $index
    """
}