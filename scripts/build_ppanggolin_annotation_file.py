import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="Build annotation file for ppanggolin")
    parser.add_argument(
        "-i", "--input", nargs="*", required=True, help="List of input GFF files"
    )
    parser.add_argument("-o", "--output", required=True, help="Output file name")

    args = parser.parse_args()

    # Process each input file
    with open(args.output, "w") as out_file:
        for file_path in args.input:
            file_name = os.path.basename(file_path)
            base_name = os.path.splitext(file_name)[0]
            out_file.write(f"{base_name}\t{file_path}\n")


if __name__ == "__main__":
    main()
