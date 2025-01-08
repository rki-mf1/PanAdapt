include {build_liftoff_db} from './modules/build_liftoff_db.nf'
include {liftoff} from './modules/liftoff.nf'
include {download_bakta_db} from './modules/download_bakta_db.nf'
include {bakta} from './modules/bakta.nf'
include {ppanggolin} from './modules/ppanggolin.nf'
include {reformat_matrix} from './modules/reformat_matrix.nf'
include {remove_stop_codons} from './modules/remove_stop_codons.nf'
include {build_gene_families} from './modules/build_gene_families.nf'
include {filter_ambiguous_sequences} from './modules/filter_ambiguous_sequences.nf'
include {remove_duplicate_sequences} from './modules/remove_duplicate_sequences.nf'
include {merge_dups_filter} from './modules/merge_dups_filter.nf'
include {filter_min_seq_count} from './modules/filter_min_seq_count.nf'
include {mafft} from './modules/mafft.nf'
include {pal2nal} from './modules/pal2nal.nf'
include {mask_ambiguities} from './modules/mask_ambiguities.nf'
include {split_reference_from_msa} from './modules/split_reference_from_msa.nf'
include {fasttree} from './modules/fasttree.nf'
include {fubar} from './modules/fubar.nf'
include {extract_fubar} from './modules/extract_fubar.nf'
include {merge_fubar_stats} from './modules/merge_fubar_stats.nf'
include {merge_gene_stats} from './modules/merge_gene_stats.nf'
include {build_site_stats} from './modules/build_site_stats.nf'
include {merge_reference} from './modules/merge_reference.nf'

workflow {

    sample_genomes = Channel.fromPath("$params.input_dir/*.fasta")
    reference_genome = Channel.value(params.reference_fasta)
    reference_gff = Channel.value(params.reference_gff)
    bakta_db = Channel.value(params.bakta_db)

    if (params.reference_fasta && params.reference_gff) {
        liftoff_db = build_liftoff_db(reference_gff)
        annotated_genomes = liftoff(liftoff_db, reference_genome, sample_genomes.mix(reference_genome))
    }
    else {
        if (params.bakta_db) {
            bakta_db = Channel.value(params.bakta_db)
        } else {
            bakta_db = download_bakta_db()
        }
        annotated_genomes = bakta(sample_genomes, bakta_db)
    }
    (ppanggolin_matrix, all_genes_nuc) = ppanggolin(annotated_genomes.collect())
    reformatted_matrix = reformat_matrix(ppanggolin_matrix)
    (gene_families_nuc, gene_families_prot, ambig_filter_counts, stop_filter_counts)  = build_gene_families(reformatted_matrix, all_genes_nuc)
    gene_families_msa = mafft(gene_families_prot.flatten())
    gene_families_nuc
        .flatten()
        .map{file->tuple(file.simpleName, file)}
        .set{gene_families_nuc_indexed}
    gene_families_msa
        .map{file->tuple(file.simpleName, file)}
        .set{gene_families_msa_indexed}

    codon_aware_msa = pal2nal(gene_families_msa_indexed.join(gene_families_nuc_indexed))
    if (params.reference_fasta && params.reference_gff) {
        (codon_aware_msa, aligned_reference) = split_reference_from_msa(codon_aware_msa)
    }
    masked_msa = mask_ambiguities(codon_aware_msa)  
    (msa_without_duplicates, dups_filter_counts) = remove_duplicate_sequences(masked_msa)
    merged_dups_filter = merge_dups_filter(dups_filter_counts.collect())
    newick_tree = fasttree(msa_without_duplicates)
    fubar_json = fubar(msa_without_duplicates.join(newick_tree))
    (fubar_per_site, fubar_stats) = extract_fubar(fubar_json)
    merged_fubar_stats = merge_fubar_stats(fubar_stats.collect())
    merged_gene_stats = merge_gene_stats(reformatted_matrix, ambig_filter_counts, stop_filter_counts, merged_dups_filter, merged_fubar_stats)
    site_stats = build_site_stats(masked_msa.join(fubar_per_site))
    if (params.reference_fasta && params.reference_gff) {
        site_stats = merge_reference(site_stats.join(aligned_reference))
    }
}   
