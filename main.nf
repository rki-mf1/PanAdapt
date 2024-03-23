input_genomes = Channel.fromPath("${params.input_genomes}/*")
reference_fasta = Channel.value(params.reference_fasta)
reference_gff = Channel.value(params.reference_gff)

params.publish_path = "results/${file(params.input_genomes).getBaseName()}"

include {seqkit_split} from './modules/seqkit_split.nf'
include {liftoff_reference} from './modules/liftoff_reference.nf'
include {liftoff} from './modules/liftoff.nf'
include {filter_invalid_orfs} from './modules/filter_invalid_orfs.nf'
include {ppanggolin} from './modules/ppanggolin.nf'
include {extract_ppanggolin_results} from './modules/extract_ppanggolin_results.nf'
include {remove_stop_codons} from './modules/remove_stop_codons.nf'
include {build_gene_families} from './modules/build_gene_families.nf'
include {sort_gene_families} from './modules/sort_gene_families.nf'
include {remove_duplicate_sequences} from './modules/remove_duplicate_sequences.nf'
include {translate_gene_families} from './modules/translate_gene_families.nf'
include {mafft} from './modules/mafft.nf'
include {pal2nal} from './modules/pal2nal.nf'
include {split_reference_from_msa} from './modules/split_reference_from_msa.nf'
include {fasttree} from './modules/fasttree.nf'
include {parnas} from './modules/parnas.nf'
include {filter_reduced_tree} from './modules/filter_reduced_tree.nf'
include {reduce_msa} from './modules/reduce_msa.nf'
include {codeml} from './modules/codeml.nf'
include {extract_codeml_parameters} from './modules/extract_codeml_parameters.nf'
include {compare_codeml_models} from './modules/compare_codeml_models.nf'
include {combine_codeml_comparisons} from './modules/combine_codeml_comparisons.nf'
include {extract_codeml_beb} from './modules/extract_codeml_beb.nf'
include {build_site_table} from './modules/build_site_table.nf'
include {build_site_beb_table} from './modules/build_site_beb_table.nf'
include {fubar} from './modules/fubar.nf'
include {build_site_fubar_table} from './modules/build_site_fubar_table.nf'
include {mask_ambiguities} from './modules/mask_ambiguities.nf'
include {filter_min_seq_count} from './modules/filter_min_seq_count.nf'
include {remove_duplicate_sequences_2} from './modules/remove_duplicate_sequences_2.nf'
include {restore_duplicate_sequences} from './modules/restore_duplicate_sequences.nf'
include {restore_duplicate_sequences_2} from './modules/restore_duplicate_sequences_2.nf'
include {filter_min_seq_count_2} from './modules/filter_min_seq_count_2.nf'
include {calculate_shannon_entropy} from './modules/calculate_shannon_entropy.nf'
include {filter_ambiguous_sequences} from './modules/filter_ambiguous_sequences.nf'
include {bakta} from './modules/bakta.nf'
include {build_refless_site_table} from './modules/build_refless_site_table.nf'
include {busted} from './modules/busted.nf'

// workflow {
//     annotated_genomes = bakta(input_genomes)
//     (matrix_csv, all_genes_fasta) = ppanggolin(annotated_genomes.collect())
//     ppanggolin_results = extract_ppanggolin_results(matrix_csv)
//     all_genes_no_stops = remove_stop_codons(all_genes_fasta)
//     gene_families = build_gene_families(matrix_csv, all_genes_no_stops)
//     gene_families_filtered = filter_ambiguous_sequences(gene_families.flatten())
//     (gene_families_without_duplicates, duplicate_map) = remove_duplicate_sequences(gene_families_filtered)
//     gene_families_min_count = filter_min_seq_count(gene_families_without_duplicates)
//     gene_families_min_count
//         .filter{file -> file.name.endsWith('.filtered')}
//         .set{gene_families_min_count}
//     gene_families_translated = translate_gene_families(gene_families_min_count)
//     gene_families_msa = mafft(gene_families_translated)
//     gene_families_min_count
//         .map{file->tuple(file.simpleName, file)}
//         .set{gene_families_min_count_indexed}
//     gene_families_msa
//         .map{file->tuple(file.simpleName, file)}
//         .set{gene_families_msa_indexed}
//     pal2nal_input = gene_families_msa_indexed.join(gene_families_min_count_indexed)
//     gene_families_codon_aware_msa = pal2nal(pal2nal_input)
//     masked_msa = mask_ambiguities(gene_families_codon_aware_msa)
//     (msa_without_duplicates, msa_duplicate_map) = remove_duplicate_sequences_2(masked_msa)
//     msa_min_count = filter_min_seq_count_2(msa_without_duplicates)
//     msa_min_count
//         .filter{index,file -> file.name.endsWith('.filtered')}
//         .set{msa_min_count}    
//     newick_tree = fasttree(msa_min_count)
//     fubar_json = fubar(msa_min_count.join(newick_tree))    
//     duplicate_map
//         .map{file->tuple(file.simpleName, file)}
//         .set{duplicate_map_indexed}
//     msa_duplicates_restored = restore_duplicate_sequences(masked_msa.join(duplicate_map_indexed))
//     site_table = build_refless_site_table(msa_duplicates_restored)
//     site_fubar_table = build_site_fubar_table(site_table.join(fubar_json))
//     site_shannon_table = calculate_shannon_entropy(site_fubar_table)

//     (matrix_csv, all_genes_fasta) = ppanggolin_v2(annotated_genomes.collect())
//     busted_json = busted(msa_min_count.join(newick_tree))   
// }


workflow {
    (annotated_reference, reference_db) = liftoff_reference(reference_fasta, reference_gff)
    annotated_genomes = liftoff(input_genomes, reference_fasta, reference_db)
    combined_annotations = annotated_genomes.mix(annotated_reference)
    filtered_annotations = filter_invalid_orfs(combined_annotations)
    (matrix_csv, all_genes_fasta) = ppanggolin(filtered_annotations.collect())
    ppanggolin_results = extract_ppanggolin_results(matrix_csv)
    all_genes_no_stops = remove_stop_codons(all_genes_fasta)
    gene_families = build_gene_families(matrix_csv, all_genes_no_stops)
    gene_families_filtered = filter_ambiguous_sequences(gene_families.flatten())
    // gene_families_sorted = sort_gene_families(gene_families.flatten())
    (gene_families_without_duplicates, duplicate_map) = remove_duplicate_sequences(gene_families_filtered)
    gene_families_min_count = filter_min_seq_count(gene_families_without_duplicates)
    gene_families_min_count
        .filter{file -> file.name.endsWith('.filtered')}
        .set{gene_families_min_count}

    gene_families_translated = translate_gene_families(gene_families_min_count)
    gene_families_msa = mafft(gene_families_translated)
    gene_families_min_count
        .map{file->tuple(file.simpleName, file)}
        .set{gene_families_min_count_indexed}
    gene_families_msa
        .map{file->tuple(file.simpleName, file)}
        .set{gene_families_msa_indexed}
    pal2nal_input = gene_families_msa_indexed.join(gene_families_min_count_indexed)
    gene_families_codon_aware_msa = pal2nal(pal2nal_input)
    masked_msa = mask_ambiguities(gene_families_codon_aware_msa)
    (msa_no_reference, aligned_reference) = split_reference_from_msa(masked_msa)
    (msa_without_duplicates, msa_duplicate_map) = remove_duplicate_sequences_2(msa_no_reference)
    msa_min_count = filter_min_seq_count_2(msa_without_duplicates)
    msa_min_count
        .filter{index,file -> file.name.endsWith('.filtered')}
        .set{msa_min_count}    
    newick_tree = fasttree(msa_min_count)
    fubar_json = fubar(msa_min_count.join(newick_tree))    
    duplicate_map
        .map{file->tuple(file.simpleName, file)}
        .set{duplicate_map_indexed}
    msa_duplicates_restored = restore_duplicate_sequences(msa_no_reference.join(duplicate_map_indexed))
    site_table = build_site_table(msa_duplicates_restored.join(aligned_reference))
    site_fubar_table = build_site_fubar_table(site_table.join(fubar_json))
    site_shannon_table = calculate_shannon_entropy(site_fubar_table)
}




//     // codeml_input = reduced_msa.join(reduced_tree_filtered).combine(params.codeml_models)
//     // (codeml_txt, codeml_rst) = codeml(codeml_input)
//     // codeml_parameters = extract_codeml_parameters(codeml_txt)
//     // codeml_parameters
//     //     .splitCsv(sep:',')
//     //     .map { items -> 
//     //         return [items[0],
//     //         items[1].toDouble(),
//     //         items[2].toInteger(),
//     //         items[3].toInteger()]
//     //     }
//     //     .groupTuple()
//     //     .set{codeml_parameters}
//     // codeml_comparison = compare_codeml_models(codeml_parameters)
//     // codeml_comparison_file = combine_codeml_comparisons(codeml_comparison.collect())
//     // codeml_beb = extract_codeml_beb(codeml_rst)
//     // codeml_beb
//     //     .filter { index, file -> file.name.endsWith('.beb') }
//     //     .set{ codeml_beb }
//     // site_beb_table = build_site_beb_table(site_table.join(codeml_beb))

//     // reduced_tree = parnas(newick_tree)
//     // reduced_tree_filtered = filter_reduced_tree(reduced_tree)
//     // reduced_tree_filtered
//     //     .filter{index, file -> file.name.endsWith('.filtered')}
//     //     .set{reduced_tree_filtered}
//     // reduced_msa = reduce_msa(msa_no_reference.join(reduced_tree_filtered))
    
// }