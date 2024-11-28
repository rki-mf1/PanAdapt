process liftoff {
    label 'liftoff'
    publishDir "${params.publish_path}/liftoff/", mode: params.publish_dir_mode
    scratch true
    
    input:
    val reference_db
    val reference_fasta
    path genome_fasta
    

    output:
    path "${genome_fasta.baseName}.gff_polished"

    script:
    """
    liftoff -db $reference_db \
            -o ${genome_fasta.baseName}.gff \
            -cds \
            -polish \
            $genome_fasta \
            $reference_fasta

    format_annotation.py -a ${genome_fasta.baseName}.gff_polished -f $genome_fasta
    sed -i '0,/CDS/{/CDS/d;}' ${genome_fasta.baseName}.gff_polished
    """
}