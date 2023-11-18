input_genomes = Channel.fromPath(params.input_genomes)
reference_fasta = Channel.value(params.reference_fasta)
reference_gff = Channel.value(params.reference_gff)

include {seqkit_split} from './modules/seqkit_split.nf'
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
include {nexus_to_newick} from './modules/nexus_to_newick.nf'
include {reduce_msa} from './modules/reduce_msa.nf'
// include {codeml} from './modules/codeml.nf'

workflow{
    split_genomes = seqkit_split(params.input_genomes)
    combined_genomes = split_genomes.flatten().mix(reference_fasta)
    annotated_genomes = liftoff(combined_genomes, reference_fasta, reference_gff)
    filtered_annotations = filter_invalid_orfs(annotated_genomes)
    (matrix_csv, all_genes_fasta) = ppanggolin(filtered_annotations.collect())
    all_genes_no_stops = remove_stop_codons(all_genes_fasta)
    gene_families = build_gene_families(matrix_csv, all_genes_no_stops)
    gene_families_sorted = sort_gene_families(gene_families.flatten())
    (gene_families_without_duplicates, duplicate_map) = remove_duplicate_sequences(gene_families_sorted)
    gene_families_translated = translate_gene_families(gene_families_without_duplicates)
    gene_families_msa = mafft(gene_families_translated)
    gene_families_without_duplicates
        .map{file->tuple(file.simpleName, file)}
        .set{gene_families_without_duplicates_indexed}
    gene_families_msa
        .map{file->tuple(file.simpleName, file)}
        .set{gene_families_msa_indexed}
    pal2nal_input = gene_families_msa_indexed.join(gene_families_without_duplicates_indexed)
    gene_families_codon_aware_msa = pal2nal(pal2nal_input)
    (msa_no_reference, aligned_reference) = split_reference_from_msa(gene_families_codon_aware_msa)
    newick_tree = fasttree(msa_no_reference)
    reduced_tree = parnas(newick_tree)
    reduced_newick_tree = nexus_to_newick(reduced_tree)
    reduced_msa = reduce_msa(msa_no_reference.join(reduced_newick_tree))
    // codeml_input = reduced_msa.join(reduced_newick_tree).combine(params.codeml_models)
    // codeml_input.view()

    // ppanggolin_results = extract_ppanggolin_results(matrix_csv)
}