sc2_genomes = Channel.fromPath(params.sc2_genomes)
reference_sequence = Channel.value(params.reference_sequence)
reference_annotation = Channel.value(params.reference_annotation)
downsample_size = Channel.value(params.downsample_size)
codeml_template = Channel.value(params.codeml_template)
codeml_models = Channel.from(params.codeml_models)

sc2_genomes_test = Channel.value(params.sc2_genomes)

process split_cds_1 {
    input:
    path genome

    output:
    path "output/*"

    script:
    """
    mkdir output
    seqkit split -i --by-id-prefix '' --out-dir output $genome
    """
}

process liftoff {
    input:
    path genome
    val reference_annotation
    val reference_sequence

    output:
    tuple path("${genome.baseName}.gff"), path(genome)

    script:
    """
    liftoff -db $reference_annotation -cds -o ${genome.baseName}.gff $genome $reference_sequence
    
    if ${params.debug}; then
        mkdir -p $baseDir/debug/liftoff
        cp ${genome.baseName}.gff $baseDir/debug/liftoff
    fi
    """
}


process remove_nonsense_sequences{
    input:
    tuple path(annotation), path(genome)

    output:
    path "output/${annotation.name}"

    script:
    """
    mkdir output
    python $baseDir/scripts/remove_nonsense_sequences.py -i $annotation -o ${annotation.baseName}.tmp
    cat ${annotation.baseName}.tmp <(echo "##FASTA") $genome > output/${annotation.name}

    if ${params.debug}; then
        mkdir -p $baseDir/debug/remove_nonsense_sequences
        cp output/${annotation.name} $baseDir/debug/remove_nonsense_sequences
    fi
    """
}

process annotate_genomes {
    input:
    path fasta_file
    val reference

    output:
    path '**.ffn'
    path '**.gff'

    script:
    """
    prokka --kingdom Viruses --proteins $reference --prefix ${fasta_file.getBaseName()} $fasta_file
    """
}

process remove_stop_codons {
    input:
    path ffn_file

    output:
    path "output/${ffn_file.name}"

    script:
    """
    mkdir output
    python $baseDir/scripts/remove_stop_codons.py -i $ffn_file -o output/${ffn_file.name}
    """
}

process split_cds_2 {
    input:
    path no_stops_file

    output:
    path "output/*"

    script:
    """
    mkdir output
    seqkit split -i --by-id-prefix '' --out-dir output $no_stops_file
    """
}

process build_pangenome {
    input:
    path gff_files

    output:
    path "pangenome/output/matrix.csv"
    path "pangenome/output/projection/*.tsv"
    path "pangenome/fasta/all_genes.fna"

    script:
    """    
    python $baseDir/scripts/build_ppanggolin_annotation_file.py -i $gff_files -o annotation.tmp
    ppanggolin annotate --cpu 12 --anno annotation.tmp --output pangenome
    ppanggolin cluster --cpu 12 -p pangenome/pangenome.h5
    ppanggolin graph --cpu 12 -p pangenome/pangenome.h5
    ppanggolin partition --cpu 12 -K 2 -p pangenome/pangenome.h5
    ppanggolin write --cpu 12 -p pangenome/pangenome.h5 --csv --projection --output pangenome/output
    ppanggolin fasta -p pangenome/pangenome.h5 --output pangenome/fasta/ --genes all
    """
}

process extract_pangenome_data {
    input:
    path pangenome_matrix

    output:
    path "results_pangenome.tsv"

    script:
    """
    awk 'BEGIN { FPAT = "([^,]+)|(\\"[^\\"]+\\")"; OFS="\\t" } { \$1=\$1; print \$0 }' $pangenome_matrix \
    | tr -d '"' | cut -f1-5 > results_pangenome.tsv
    """
}

process build_families {
    input:
    path pan_matrix
    path split_cds_files

    output:
    path "output/*.fasta"
    path "output/*.ref"

    script:
    """
    mkdir output
    cat $split_cds_files > combined_cds.tmp
    python $baseDir/scripts/build_families.py -m $pan_matrix -s combined_cds.tmp -r ${params.ref_id} -o output

    if ${params.debug}; then
        mkdir -p $baseDir/debug/build_families
        cp output/*.fasta $baseDir/debug/build_families
        cp output/*.ref $baseDir/debug/build_families
    fi
    """
}

process sort_families {
    input:
    tuple val(index), path(family)

    output:
    tuple val(index), path("output/${family.name}")

    script:
    """
    mkdir output
    seqkit sort $family > output/${family.name}
    
    mkdir -p $baseDir/debug/sorted_families
    cp output/${family.name} $baseDir/debug/sorted_families
    """
}

process filter_sequences_by_length {
    input:
    tuple val(index), path(family)

    output:
    tuple val(index), path("output/${family.baseName}.fasta")

    script:
    """
    mkdir output
    python $baseDir/scripts/filter_sequences_by_length.py -i $family -o output/${family.baseName}.fasta
    """
}

process remove_duplicate_sequences {
    input:
    tuple val(index), path(filtered_family)

    output:
    tuple val(index), path("output/${filtered_family.name}")
    tuple val(index), path("output/${filtered_family.getBaseName()}.json")

    script:
    """
    mkdir output
    python $baseDir/scripts/remove_duplicate_sequences.py -i $filtered_family -o output/${filtered_family.name} -j output/${filtered_family.getBaseName()}.json
    """
}


process count_filtered_sequences {
    input:
    tuple val(index), path(no_dups_family), path(ref_seq)

    output:
    tuple val(index), path("output/${no_dups_family.baseName}.*")
    path ".command.out" 

    script:
    """
    mkdir output
    count=`grep -c '^>' $no_dups_family`
    echo -e "${no_dups_family.baseName}\t\$count"
    if [ \$count -ge 2 ]; then 
        cp $no_dups_family output/${no_dups_family.name}
        cat $ref_seq >> output/${no_dups_family.name}
    else 
        touch output/${no_dups_family.baseName}.filter
    fi
    """
}

process write_filtered_counts_file{
    input:
    val filtered_counts

    output:
    path "filtered_counts.tsv"

    script:
    def filtered_counts_str = filtered_counts.join(' ')
    """
    echo -e "Gene\tNo. filtered" > filtered_counts.tsv
    cat $filtered_counts_str >> filtered_counts.tsv
    """
}


process translate_families {
    input:
    tuple val(index), path(family)

    output:
    tuple val(index), path("output/${family.name}")

    script:
    """
    mkdir output
    transeq -sequence $family -outseq output/${family.name}
    sed -i 's/_1\$//' output/${family.name}
    """
}

process build_msa{
    input:
    tuple val(index), path(family)

    output:
    tuple val(index), path("output/${family.baseName}.msa")

    script:
    """
    mkdir output
    mafft --auto --preservecase $family > output/${family.baseName}.msa
    """
}

process build_codon_aware_msa{
    input:
    tuple val(index), path(msa), path(family)

    output:
    tuple val(index), path("output/${family.name}")

    script:
    """
    mkdir output
    pal2nal.pl -output fasta $msa $family > output/${family.name}

    if ${params.debug}; then
    mkdir -p $baseDir/debug/build_codon_aware_msa
    cp output/${family.name} $baseDir/debug/build_codon_aware_msa
    fi
    """
}

process mask_ambiguities_with_consensus{
    input:
    tuple val(index), path(msa)

    output:
    tuple val(index), path("output/${msa.name}")

    script:
    """
    mkdir output
    python $baseDir/scripts/test.py -i $msa -o output/${msa.name}

    if ${params.debug}; then
    mkdir -p $baseDir/debug/mask_ambiguities_with_consensus
    cp output/${msa.name} $baseDir/debug/mask_ambiguities_with_consensus
    fi
    """
}

process split_ref_from_msa{
    input:
    tuple val(index), path(msa)

    output:
    tuple val(index), path("output/${msa.name}")
    tuple val(index), path("output/${msa.baseName}.ref")

    """
    mkdir output
    python $baseDir/scripts/split_ref_from_msa.py -i $msa -r ${params.ref_id} -o output/${msa.name} -u output/${msa.baseName}.ref

    if ${params.debug}; then
        mkdir -p $baseDir/debug/split_ref_from_msa
        cp output/${msa.name} $baseDir/debug/split_ref_from_msa
        cp output/${msa.baseName}.ref $baseDir/debug/split_ref_from_msa
    fi
    """
}

process build_newick_tree {
    input:
    tuple val(index), path(msa)

    output:
    tuple val(index), path("${msa.baseName}.nwk")

    script:
    """
    fasttree -nt $msa > ${msa.baseName}.nwk
    """
}

process downsample_newick_trees {
    input:
    tuple val(index), path(newick_tree)
    val downsample_size

    output:
    tuple val(index), path("output/${newick_tree.name}")

    script:
    """
    mkdir output
    default_subtree_size=25
    seq_count=`grep -o ',' $newick_tree | wc -l`
    if [ "\$seq_count" -le "$downsample_size" ]; then
        cp $newick_tree output/${newick_tree.name}
    else
        parnas -t $newick_tree --cover --radius 0.001 --subtree output/${newick_tree.name}
        sed -i "s/'//g" output/${newick_tree.name}
    fi
    """
}

process reduce_msa {
    input:
    tuple val(index), path(msa), path(newick_tree)

    output:
    tuple val (index), path("output/${msa.name}"), path("output/${newick_tree.name}")

    script:
    """
    mkdir output
    if head -n 1 $newick_tree | grep -q '^#NEXUS'; then
        python $baseDir/scripts/reduce_msa.py -i $msa -t $newick_tree -o output/${msa.name} -r output/${newick_tree.name}
    else
        cp $msa output/${msa.name}
        cp $newick_tree output/${newick_tree.name}
    fi

    if ${params.debug}; then
        mkdir -p $baseDir/debug/reduce_msa
        cp output/${msa.name} $baseDir/debug/reduce_msa
    fi
    """
}

process run_codeml_model{
    input:
    tuple val(index), path(msa), path(newick_tree), val(codeml_model)
    val codeml_template

    output:
    tuple val(index), path("${msa.baseName}.txt"), path("rst"), val(codeml_model)

    script:
    """
    sed 's|stewart\\.aa|${msa}|g;
    s|stewart\\.trees|${newick_tree}|g;
    s|out_M0\\.txt|${msa.baseName}.txt|g;
    s|NSsites_placeholder|${codeml_model}|g;' $codeml_template > codeml.ctl
    codeml codeml.ctl
    """
}

process extract_codeml_parameters{
    input:
    tuple val(index), path(codeml_output), path(rst), val(codeml_model)

    output:
    stdout

    script:
    """
    lnL=`grep 'lnL' $codeml_output | sed -n 's/.*\\: \\{0,\\}\\([-0-9.]\\+\\).*/\\1/p'`
    np=`grep 'lnL' $codeml_output | sed -n 's/.*np: \\{0,\\}\\([0-9]\\+\\).*/\\1/p'`
    echo "$index,\$lnL,\$np,$codeml_model"
    """
}

process compare_model_likelihood{
    input:
    tuple val(index), val(lnL), val(np), val(codeml_model), path(family)

    output:
    path ".command.out"

    script:
    def lnL_str = lnL.join(' ')
    def np_str = np.join(' ')
    def codeml_model_str = codeml_model.join(' ')

    """
    python $baseDir/scripts/compare_model_likelihood2.py -l $lnL_str -p $np_str -m $codeml_model_str -g ${family.baseName}
    """
}

process write_comparison_file{
    input:
    val model_comparisons

    output:
    path "model_comparison.tsv"

    script:
    def model_comparisons_str = model_comparisons.join(' ')
    """
    echo -e "Gene\tMO_vs_M1_p_value\tM1_vs_M2_p_value\tM7_vs_M8_p_value\tPreferred_Basic_Model\tPreferred_Site_Model" > model_comparison.tsv
    cat $model_comparisons_str >> model_comparison.tsv
    """
}

process merge_results {
    input:
    val pan_results
    val filtered_counts
    val model_comparison_file

    output:
    path "merged_results.tsv"

    script:
    """
    python $baseDir/scripts/merge_results.py $pan_results $filtered_counts $model_comparison_file -o merged_results.tsv

    if ${params.debug}; then
        mkdir -p $baseDir/debug/merge_results
        cp merged_results.tsv $baseDir/debug/merge_results
    fi
    """
}

process extract_beb_table {
    input:
    tuple val(index), path(codeml_output), path(rst), val(codeml_model)

    output:
    tuple val(index), path("output/${codeml_output.baseName}.*")

    script:
    """
    mkdir output
    if [ "$codeml_model" == '8' ]; then 
        python $baseDir/scripts/extract_beb_table.py -i $rst -o output/${codeml_output.baseName}.beb
    else 
        touch output/${codeml_output.baseName}.filter
    fi

    if ${params.debug}; then
        mkdir -p $baseDir/debug/extract_beb_table
        cp output/${codeml_output.baseName}.* $baseDir/debug/extract_beb_table
    fi
    """
}

process add_reference_to_msa {
    input:
    tuple val(index), path(msa), path(tree), path(ref_seq)

    output:
    tuple val(index), path("output/${msa.name}")

    script:
    """
    mkdir output
    sed 's/"-"/""/g' $msa > 1.tmp
    cat $ref_seq >> 1.tmp
    transeq -sequence 1.tmp -outseq 2.tmp
    sed -i 's/_1\$//' 2.tmp
    mafft --auto --preservecase 2.tmp > 3.tmp
    pal2nal.pl -output fasta 3.tmp 1.tmp > output/${msa.name}

    if ${params.debug}; then
        mkdir -p $baseDir/debug/add_reference_to_msa
        cp output/${msa.name} $baseDir/debug/add_reference_to_msa
    fi
    """
}

process initialize_per_site_table {
    input:
    tuple val(index), path(msa), path(tree), path(ref_seq)

    output:
    tuple val(index), path("${msa.baseName}.tsv")

    script:
    """
    echo 1
    cat $msa $ref_seq > ${msa.baseName}.tmp
    python $baseDir/scripts/initialize_per_site_table.py -i ${msa.baseName}.tmp -r ${params.ref_id} -o ${msa.baseName}.tsv
    
    if ${params.debug}; then
        mkdir -p $baseDir/debug/initialize_per_site_table
        cp ${msa.baseName}.tsv $baseDir/debug/initialize_per_site_table
    fi
    """
}

process join_per_site_and_beb {
    input:
    tuple val(index), path(per_site_table), path(beb_table)

    output:
    tuple val(index), path("output/${per_site_table.name}")

    script:
    """
    mkdir output 
    python $baseDir/scripts/join_per_site_and_beb.py -i $per_site_table -b $beb_table -o output/${per_site_table.name}
    
    if ${params.debug}; then
        mkdir -p $baseDir/debug/join_per_site_and_beb
        cp output/${per_site_table.name} $baseDir/debug/join_per_site_and_beb
    fi
    """
}

process visualize_beb {
    input:
    tuple val(index), path(per_site_table)

    output:
    tuple val(index), path("${per_site_table.baseName}.png")

    script:
    """
    echo 15
    python $baseDir/scripts/visualize_beb.py -i $per_site_table -o ${per_site_table.baseName}.png

    if ${params.debug}; then
        mkdir -p $baseDir/debug/visualize_beb
        cp ${per_site_table.baseName}.png $baseDir/debug/visualize_beb
    fi
    """
}

process extract_positively_selected_sites {
    input:
    tuple val(index), path(per_site_table)
    path sc2_genomes

    output:
    tuple val(index), path("output/${per_site_table.name}")

    script:
    """
    mkdir output
    python $baseDir/scripts/extract_positively_selected_sites.py -i $per_site_table -o output/${per_site_table.name}
    mkdir -p $baseDir/ps_sites/${sc2_genomes.baseName}
    cp output/${per_site_table.name} $baseDir/ps_sites/${sc2_genomes.baseName}
    
    if ${params.debug}; then
        mkdir -p $baseDir/debug/extract_positively_selected_sites
        cp output/${per_site_table.name} $baseDir/debug/extract_positively_selected_sites
    fi
    """
}

workflow {
    random_sample = split_cds_1(sc2_genomes)
    random_sample
        .flatten()
        .mix(reference_sequence)
        .set{random_sample}
    gffs_and_genomes = liftoff(random_sample.flatten(), reference_annotation, reference_sequence)
    nonsense_removed = remove_nonsense_sequences(gffs_and_genomes)
   (pan_matrix, pan_projections, ffn_files) = build_pangenome(nonsense_removed.collect())
    pan_results = extract_pangenome_data(pan_matrix)
    ffn_stopless_files = remove_stop_codons(ffn_files)
    ffn_split_files = split_cds_2(ffn_stopless_files)
    (families, ref_seqs) = build_families(pan_matrix, ffn_split_files.collect())
    counter = { int i=0; return { item -> [++i, item] } }()
    families = families.flatten().map(counter)
    counter2 = { int i=0; return { item -> [++i, item] } }()
    ref_seqs = ref_seqs.flatten().map(counter2)
    sorted_families = sort_families(families)
    filtered_families = filter_sequences_by_length(sorted_families)    
    (no_dups_families, dups_jsons) = remove_duplicate_sequences(filtered_families)
    families_with_refs = no_dups_families.join(ref_seqs)
    (min_count_families, seq_counts) = count_filtered_sequences(families_with_refs)
    filtered_counts = write_filtered_counts_file(seq_counts.collect())
    min_count_families
        .filter { index, file -> file.name.endsWith('.fasta') }
        .set{ min_count_families }
    translated_families = translate_families(min_count_families)
    msa_families = build_msa(translated_families)
    joined_pla2nal = msa_families.join(min_count_families)
    msa_codon_aware_families = build_codon_aware_msa(joined_pla2nal)
    msa_codon_aware_families = mask_ambiguities_with_consensus(msa_codon_aware_families)
    (refless_msas, aligned_refs) = split_ref_from_msa(msa_codon_aware_families)
    newick_trees = build_newick_tree(refless_msas)
    downsampled_trees = downsample_newick_trees(newick_trees, downsample_size)
    msas_and_trees = refless_msas.join(downsampled_trees)
    reduced_msas_and_trees = reduce_msa(msas_and_trees)
    reduced_msas_and_trees.join(aligned_refs).set{msas_and_aligned_refs}    
    per_site_tables = initialize_per_site_table(msas_and_aligned_refs)
    codeml_inputs = reduced_msas_and_trees.combine(codeml_models)
    codeml_results = run_codeml_model(codeml_inputs, codeml_template)
    codeml_parameters = extract_codeml_parameters(codeml_results)
    codeml_parameters
        .splitCsv(sep:',')
        .map { items -> 
            return [items[0].toInteger(),
            items[1].toDouble(),
            items[2].toInteger(),
            items[3].toInteger()]
        }
        .groupTuple()
        .set{ codeml_parameters }
    codeml_parameters = codeml_parameters.join(msa_codon_aware_families)
    model_comparisons = compare_model_likelihood(codeml_parameters)
    model_comparison_file = write_comparison_file(model_comparisons.collect())
    merged_results = merge_results(pan_results, filtered_counts, model_comparison_file)
    beb_tables = extract_beb_table(codeml_results)
    beb_tables
        .filter { index, file -> file.name.endsWith('.beb') }
        .set{ beb_tables }
    beb_tables.view()
    per_site_tables = per_site_tables.join(beb_tables)
    per_site_tables = join_per_site_and_beb(per_site_tables)
    beb_visualizations = visualize_beb(per_site_tables)
    ps_site_tables = extract_positively_selected_sites(per_site_tables, sc2_genomes_test)
}