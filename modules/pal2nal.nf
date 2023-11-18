process pal2nal {
    publishDir "${params.output}/pal2nal", mode: params.publish_dir_mode

    input:
    tuple val(index), path(gene_family_msa), path(gene_family)

    output:
    tuple val(index), path("${gene_family_msa.simpleName}.pal2nal")

    script:
    """
    pal2nal.pl -output fasta $gene_family_msa $gene_family > ${gene_family_msa.simpleName}.pal2nal
    """
}