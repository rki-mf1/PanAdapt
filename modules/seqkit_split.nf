process seqkit_split {
    label 'seqkit'
    publishDir "${params.publish_path}/seqkit_split/${fasta.baseName}", mode: params.publish_dir_mode
    
    input:
    path fasta

    output:
    path "**.fasta"

    script:
    """
    seqkit split -j $task.cpus -i --by-id-prefix '' $fasta
    """
}