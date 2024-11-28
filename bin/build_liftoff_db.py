#!/usr/bin/env python
import argparse
import gffutils


def build_liftoff_db(gff_file, db):
    gffutils.create_db(gff_file, db, merge_strategy="create_unique", force=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    build_liftoff_db(args.input, args.output)


if __name__ == "__main__":
    main()
