process liftoff_reference {
    label 'liftoff'
    publishDir "${params.publish_path}/liftoff_reference/", mode: params.publish_dir_mode
    
    input:
    path reference_fasta
    path reference_gff

    output:
    path "${reference_fasta.baseName}.gff"
    path "${reference_gff}_db"

    script:
    """
    liftoff -g $reference_gff \
            -cds \
            -o ${reference_fasta.baseName}.tmp \
            $reference_fasta \
            $reference_fasta
    cat ${reference_fasta.baseName}.tmp <(echo "##FASTA") $reference_fasta > ${reference_fasta.baseName}.gff
    """
}