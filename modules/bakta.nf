process bakta {
    label 'bakta'
    publishDir "${params.publish_path}/bakta/${genome_fasta.baseName}", mode: params.publish_dir_mode

    input:
    path genome_fasta

    output:
    path "${genome_fasta.baseName}.gff3"

    script:
    """
    bakta \
    --db ${params.bakta_db} \
    --threads $task.cpus \
    --skip-trna \
    --skip-tmrna \
    --skip-rrna \
    --skip-ncrna \
    --skip-ncrna-region \
    --skip-crispr \
    --skip-pseudo \
    --skip-gap \
    --skip-ori \
    --skip-plot \
    $genome_fasta
    """
}