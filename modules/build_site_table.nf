process build_site_table {
    label 'generic_small'
    publishDir "${params.publish_path}/build_site_table", mode: params.publish_dir_mode
    scratch true
    
    input:
    tuple val(index), path(msa), path(ref_seq)

    output:
    tuple val(index), path("${index}.tsv")

    script:
    """
    cat $msa $ref_seq > ${msa.baseName}.tmp
    build_site_table.py -i ${msa.baseName}.tmp -r ${params.ref_id} -o ${index}.tsv
    """
}