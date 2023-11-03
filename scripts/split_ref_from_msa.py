import argparse
from utils import read_fasta, write_fasta

def split_ref_from_msa(input_file, ref_id, output_file, ref_output_file):
    with open(output_file, 'w') as outfile, \
         open(ref_output_file, 'w') as ref_outfile:
        for identifier, sequence in read_fasta(input_file):
            if ref_id in identifier:
                write_fasta(identifier, sequence, ref_outfile)
            else:
                write_fasta(identifier, sequence, outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split reference sequence from MSA based on ref_id")
    parser.add_argument("-i", "--input", required=True, help="Path to the input MSA file")
    parser.add_argument("-r", "--ref_id", required=True, help="Reference ID to search for in the MSA headers")
    parser.add_argument("-o", "--output", required=True, help="Path to the output MSA file without reference sequence")
    parser.add_argument("-u", "--ref_output", required=True, help="Path to the output file containing the reference sequence")

    args = parser.parse_args()

    split_ref_from_msa(args.input, args.ref_id, args.output, args.ref_output)
