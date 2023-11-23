#!/usr/bin/env python
import argparse
import re

def extract_parameters(codeml_output_file):
    lnL, np = None, None
    with open(codeml_output_file, 'r') as file:
        for line in file:
            if 'lnL' in line:
                lnL_match = re.search(r':\s*(-?\d+\.\d+)', line)
                np_match = re.search(r'np:\s*(\d+)', line)
                if lnL_match:
                    lnL = lnL_match.group(1)
                if np_match:
                    np = np_match.group(1)
                break  # Assuming the first occurrence contains the needed data

    return lnL, np

def main():
    parser = argparse.ArgumentParser(description='Extract parameters from codeml output.')
    parser.add_argument('codeml_output', type=str, help='Path to the codeml output file.')
    parser.add_argument('index', type=str, help='Index value.')
    parser.add_argument('codeml_model', type=str, help='Codeml model.')

    args = parser.parse_args()

    lnL, np = extract_parameters(args.codeml_output)
    if lnL is not None and np is not None:
        print(f"{args.index},{lnL},{np},{args.codeml_model}")
    else:
        print(f"Error extracting parameters from {args.codeml_output}")

if __name__ == "__main__":
    main()
