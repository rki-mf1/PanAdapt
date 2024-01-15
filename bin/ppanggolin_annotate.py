#!/usr/bin/env python
import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")

    args = parser.parse_args()
    with open(args.input) as infile:
        file_paths = infile.read().splitlines()

    with open(args.output, "w") as out_file:
        for file_path in file_paths:
            file_name = os.path.basename(file_path)
            base_name = os.path.splitext(file_name)[0]
            out_file.write(f"{base_name}\t{file_path}\n")


if __name__ == "__main__":
    main()
