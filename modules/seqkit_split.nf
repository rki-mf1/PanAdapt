process seqkit_split {
    publishDir "${params.publish_path}/seqkit_split/${fasta.baseName}", mode: params.publish_dir_mode
    
    input:
    path fasta

    output:
    path "**.fasta"

    script:
    """
    seqkit split -i --by-id-prefix '' $fasta
    """
}