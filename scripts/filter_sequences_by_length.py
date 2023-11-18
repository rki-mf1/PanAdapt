import argparse
from utils import read_fasta, write_fasta


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter sequences in a fasta file based on their length compared to the mean."
    )
    parser.add_argument("-i", "--infile", required=True, help="Input fasta file.")
    parser.add_argument(
        "-o",
        "--outfile",
        required=True,
        help="Output fasta file after filtering sequences.",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.1,
        help="Threshold for filtering sequences based on deviation from mean length. Default is 0.1 (10%).",
    )
    return parser.parse_args()


def filter_sequences_by_length(infile, outfile, threshold_percentage=10):
    original_sequences = list(read_fasta(infile))
    sequences = original_sequences.copy()

    initial_mean_length = (
        sum(len(seq) for _, seq in original_sequences) / len(original_sequences)
        if original_sequences
        else 0
    )

    while True:
        lengths = [len(seq) for _, seq in sequences]
        mean_length = sum(lengths) / len(lengths) if lengths else 0
        threshold = threshold_percentage / 100 * mean_length
        lower_bound = mean_length - threshold
        upper_bound = mean_length + threshold

        # Check if all sequences are within the threshold
        if all(lower_bound <= length <= upper_bound for length in lengths):
            break

        # Find the index of the most atypical sequence
        deviations = [abs(length - mean_length) for length in lengths]
        atypical_index = deviations.index(max(deviations))

        # Remove the most atypical sequence
        del sequences[atypical_index]

    with open(outfile, "w") as out_fasta:
        for identifier, seq in sequences:
            write_fasta(identifier, seq, out_fasta)

    print(f"Initial sequences: {len(original_sequences)}")
    print(f"Sequences after filtering: {len(sequences)}")
    print(f"Percentage retained: {len(sequences)/len(original_sequences)*100:.2f}%")
    print(f"Initial mean length: {initial_mean_length:.2f}")
    print(f"Final mean length: {mean_length:.2f}")


def main():
    args = parse_args()
    filter_sequences_by_length(args.infile, args.outfile, args.threshold)


if __name__ == "__main__":
    main()