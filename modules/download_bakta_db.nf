 process download_bakta_db {
    label 'bakta'
    publishDir "${params.publish_path}/download_bakta_db", mode: params.publish_dir_mode

    output:
    path "db-light"

    script:
    """
    bakta_db download --output . --type light
    """
 }