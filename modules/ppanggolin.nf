process ppanggolin {
    label 'ppanggolin'
    publishDir "${params.publish_path}/ppanggolin", mode: params.publish_dir_mode

    input:
    path gff_files

    output:
    file "pangenome/matrix.csv"
    file "all_genes/all_genes.fna"

    script:
    """
    echo $gff_files | tr ' ' '\n' > test.tmp    
    ppanggolin_annotate.py -i test.tmp -o annotation.tmp
    ppanggolin all --anno annotation.tmp --output pangenome
    ppanggolin fasta -p pangenome/pangenome.h5 --genes all --output all_genes/
    """
}
