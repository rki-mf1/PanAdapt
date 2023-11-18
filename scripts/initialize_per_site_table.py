import argparse
from utils import read_fasta

def build_per_site_table(msa_file, output_file, ref_id=None):
    identifiers = []
    sequences = []
    for identifier, seq in read_fasta(msa_file):
        identifiers.append(identifier)
        sequences.append(seq)

    # Ensure all sequences are of the same length
    if not all(len(seq) == len(sequences[0]) for seq in sequences):
        raise ValueError("All sequences in the MSA must have the same length")

    # Ensure sequence length is a multiple of 3 (codon length)
    if len(sequences[0]) % 3 != 0:
        raise ValueError("Sequence length must be a multiple of 3 for codon-based positions")

    ref_index = None
    if ref_id:
        for idx, identifier in enumerate(identifiers):
            if ref_id in identifier:
                ref_index = idx
                break
        # If ref_id is not found, use the first sequence as the reference
        if ref_index is None:
            ref_index = 0

    with open(output_file, 'w') as out:
        # Write header
        headers = ["Position", "Position_in_reference"]
        if ref_index is not None:
            headers.append("Ref_codons")
        headers.extend([f"Seq_{i+1}" for i, _ in enumerate(identifiers)])
        out.write("\t".join(headers) + "\n")
        
        # Write per-codon data
        ref_pos = 0
        for pos in range(0, len(sequences[0]), 3):
            codons = [seq[pos:pos+3] for seq in sequences]
            if ref_index is not None:
                ref_codon = codons[ref_index]
                if ref_codon != '---':
                    ref_pos += 1
                out.write(f"{(pos//3)+1}\t{ref_pos}\t{ref_codon}")
            else:
                out.write(f"{(pos//3)+1}\t{ref_pos}")
            out.write("\t" + "\t".join(codons) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build a per-site table from a multiple sequence alignment.")
    parser.add_argument('-i', '--input', required=True, help="Input MSA file.")
    parser.add_argument('-o', '--output', required=True, help="Output table file.")
    parser.add_argument('-r', '--ref_id', help="Substring of the reference ID. If present in any MSA identifier, its column will be renamed to 'Ref_codons' and moved next to 'Positions'.")
    args = parser.parse_args()

    build_per_site_table(args.input, args.output, args.ref_id)