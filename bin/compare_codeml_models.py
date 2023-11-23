#!/usr/bin/env python
import argparse
from scipy.stats import chi2


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--lnL", required=True, nargs="+", type=float)
    parser.add_argument("-p", "--np", required=True, nargs="+", type=int)
    parser.add_argument("-m", "--codeml_model", required=True, nargs="+", type=int)
    parser.add_argument("-g", "--gene", required=True, type=str)
    return parser.parse_args()

def calculate_LRT_and_dof(null_lnL, alt_lnL, null_np, alt_np):
    LRT = 2 * (alt_lnL - null_lnL)
    dof = alt_np - null_np
    return LRT, dof

def calculate_p_value(chi2_value, dof):
    return chi2.sf(chi2_value, dof)

def main():
    args = parse_args()
    lnL_dict = {k:v for k, v in zip(args.codeml_model, args.lnL)}
    np_dict = {k:v for k, v in zip(args.codeml_model, args.np)}
    diff_M0vsM1, dof_M0vsM1 = calculate_LRT_and_dof(lnL_dict[0], lnL_dict[1], np_dict[0], np_dict[1])
    p_value_M0vsM1 = calculate_p_value(diff_M0vsM1, dof_M0vsM1)

    diff_M1vsM2, dof_M1vsM2 = calculate_LRT_and_dof(lnL_dict[1], lnL_dict[2], np_dict[1], np_dict[2])
    p_value_M1vsM2 = calculate_p_value(diff_M1vsM2, dof_M1vsM2)

    diff_M7vsM8, dof_M7vsM8 = calculate_LRT_and_dof(lnL_dict[7], lnL_dict[8], np_dict[7], np_dict[8])
    p_value_M7vsM8 = calculate_p_value(diff_M7vsM8, dof_M7vsM8)

    preferred_basic_model = "M0"
    if p_value_M1vsM2 < 0.05:
        preferred_basic_model = "M2"
    elif p_value_M0vsM1 < 0.05:
        preferred_basic_model = "M1"

    preferred_site_model = "M8" if p_value_M7vsM8 < 0.05 else "M7"
    
    values = [
        args.gene,
        "{:.3f}".format(p_value_M0vsM1),
        "{:.3f}".format(p_value_M1vsM2),
        "{:.3f}".format(p_value_M7vsM8),
        preferred_basic_model,
        preferred_site_model,
    ]
    print("\t".join(values))

if __name__ == "__main__":
    main()
