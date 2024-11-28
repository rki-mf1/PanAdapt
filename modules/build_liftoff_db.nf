process build_liftoff_db {
    label 'liftoff'
    publishDir "${params.publish_path}/build_liftoff_db/", mode: params.publish_dir_mode

    
    input:
    path reference_gff
    

    output:
    path "${reference_gff.baseName}.db"

    script:
    """
    build_liftoff_db.py -i $reference_gff -o ${reference_gff.baseName}.db
    """
}