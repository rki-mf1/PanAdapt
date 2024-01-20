process ppanggolin {
    label 'generic_large'
    publishDir "${params.publish_path}/ppanggolin", mode: params.publish_dir_mode

    input:
    path gff_files

    output:
    file "matrix/matrix.csv"
    file "all_genes/all_genes.fna"

    script:
    """
    echo $gff_files | tr ' ' '\n' > test.tmp    
    ppanggolin_annotate.py -i test.tmp -o annotation.tmp
    ppanggolin annotate --cpu $task.cpus --anno annotation.tmp --output pangenome
    ppanggolin cluster --cpu $task.cpus -p pangenome/pangenome.h5
    ppanggolin graph --cpu $task.cpus -p pangenome/pangenome.h5
    ppanggolin partition --cpu $task.cpus -K 2 -p pangenome/pangenome.h5
    ppanggolin write --cpu $task.cpus -p pangenome/pangenome.h5 --csv -o matrix/
    ppanggolin fasta -p pangenome/pangenome.h5 --genes all --output all_genes/
    """
}
