process ppanggolin {
    publishDir "${params.publish_path}/ppanggolin", mode: params.publish_dir_mode

    input:
    path gff_files

    output:
    file "matrix/matrix.csv"
    file "all_genes/all_genes.fna"

    script:
    """    
    ppanggolin_annotate.py -i $gff_files -o annotation.tmp
    ppanggolin annotate --cpu 12 --anno annotation.tmp --output pangenome
    ppanggolin cluster --cpu 12 -p pangenome/pangenome.h5
    ppanggolin graph --cpu 12 -p pangenome/pangenome.h5
    ppanggolin partition --cpu 12 -K 2 -p pangenome/pangenome.h5
    ppanggolin write --cpu 12 -p pangenome/pangenome.h5 --csv -o matrix/
    ppanggolin fasta -p pangenome/pangenome.h5 --genes all --output all_genes/
    """
}
