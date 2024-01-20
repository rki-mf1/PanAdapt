process liftoff {
    label 'generic_small'
    publishDir "${params.publish_path}/liftoff/${genome_fasta.baseName}", mode: params.publish_dir_mode
    scratch true
    
    input:
    path genome_fasta
    val reference_fasta
    val reference_gff

    output:
    path "${genome_fasta.baseName}.gff"

    script:
    """
    liftoff -db $reference_gff \
            -cds \
            -o ${genome_fasta.baseName}.tmp \
            $genome_fasta \
            $reference_fasta
    cat ${genome_fasta.baseName}.tmp <(echo "##FASTA") $genome_fasta > ${genome_fasta.baseName}.gff
    """
}