import random


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


def get_random_sequences(infile, total_sequences, n):
    random_indices = random.sample(range(total_sequences), n)
    random_sequences = []
    for i, (identifier, sequence) in enumerate(read_fasta(infile)):
        if i % 1000 == 0:
            print(f"Processed {i} sequences...")
        if i in random_indices:
            random_sequences.append((identifier, sequence))
    print(f"Processed {i + 1} sequences in total.")
    return random_sequences
