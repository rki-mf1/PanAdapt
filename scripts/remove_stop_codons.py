import argparse


def read_fasta(infile):
    with open(infile, "r") as file:
        identifier = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(sequence)
                identifier = line[1:]
                sequence = []
            else:
                sequence.append(line)

        if identifier is not None:
            yield identifier, "".join(sequence)


def write_fasta(identifier, sequence, outfile, line_length=60):
    print(f">{identifier}", file=outfile)
    formatted_sequence = "\n".join(
        [sequence[i : i + line_length] for i in range(0, len(sequence), line_length)]
    )
    print(formatted_sequence, file=outfile)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    return parser.parse_args()


args = parse_args()
with open(args.outfile, "w") as outfile:
    for identifier, sequence in read_fasta(args.infile):
        if sequence[-3:] in ["TAA", "TAG", "TGA"]:
            sequence = sequence[:-3]
        write_fasta(identifier, sequence, outfile)
