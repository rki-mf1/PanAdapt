process remove_stop_codons {
    label 'generic_small'
    publishDir "${params.publish_path}/remove_stop_codons", mode: params.publish_dir_mode

    input:
    path fasta

    output:
    path "${fasta.baseName}.fasta"

    script:
    """
    remove_stop_codons.py -i $fasta -o ${fasta.baseName}.fasta
    """
}