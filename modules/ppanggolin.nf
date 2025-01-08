process ppanggolin {
    label 'ppanggolin'
    publishDir "${params.publish_path}/ppanggolin", mode: params.publish_dir_mode

    input:
    path gff_files

    output:
    file "pangenome/matrix.csv"
    file "all_genes/all_genes.fna"
    file "all_prots/all_protein_genes.faa"

    script:
    """
    echo $gff_files | tr ' ' '\n' > test.tmp    
    ppanggolin_annotate.py -i test.tmp -o ppanggolin_annotations.tsv
    ppanggolin all --anno ppanggolin_annotations.tsv --output pangenome
    ppanggolin fasta -p pangenome/pangenome.h5 --genes all --output all_genes/
    ppanggolin fasta -p pangenome/pangenome.h5 --protein all --output all_prots/
    """
}
