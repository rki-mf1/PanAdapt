import argparse
import pandas as pd
from scipy.stats import chi2


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", required=True)
    parser.add_argument("-o", "--outfile", required=True)
    return parser.parse_args()


def calculate_LRT_and_dof(df, null_model, alt_model):
    null_lnL = df[df["Model"] == null_model]["lnL"].values[0]
    alt_lnL = df[df["Model"] == alt_model]["lnL"].values[0]
    print(null_lnL, alt_lnL)
    LRT = 2 * (alt_lnL - null_lnL)
    null_np = df[df["Model"] == null_model]["np"].values[0]
    alt_np = df[df["Model"] == alt_model]["np"].values[0]
    dof = alt_np - null_np
    return LRT, dof


def calculate_p_value(chi2_value, dof):
    return chi2.sf(chi2_value, dof)


def main():
    args = parse_args()
    data = pd.read_csv(args.infile, sep="\t")
    results = []
    for gene, group in data.groupby("Gene"):
        print(group)
        diff_M0vsM1, dof_M0vsM1 = calculate_LRT_and_dof(group, 0, 1)
        p_value_M0vsM1 = calculate_p_value(diff_M0vsM1, dof_M0vsM1)

        diff_M1vsM2, dof_M1vsM2 = calculate_LRT_and_dof(group, 1, 2)
        p_value_M1vsM2 = calculate_p_value(diff_M1vsM2, dof_M1vsM2)

        diff_M7vsM8, dof_M7vsM8 = calculate_LRT_and_dof(group, 7, 8)
        p_value_M7vsM8 = calculate_p_value(diff_M7vsM8, dof_M7vsM8)

        preferred_basic_model = "M0"
        if p_value_M1vsM2 < 0.05:
            preferred_basic_model = "M2"
        elif p_value_M0vsM1 < 0.05:
            preferred_basic_model = "M1"

        preferred_site_model = "M8" if p_value_M7vsM8 < 0.05 else "M7"

        results.append(
            {
                "Gene": gene,
                "M0_vs_M1_p_value": "{:.3f}".format(p_value_M0vsM1),
                "M1_vs_M2_p_value": "{:.3f}".format(p_value_M1vsM2),
                "M7_vs_M8_p_value": "{:.3f}".format(p_value_M7vsM8),
                "Preferred_Basic_Model": preferred_basic_model,
                "Preferred_Site_Model": preferred_site_model,
            }
        )

    results_df = pd.DataFrame(results)
    results_df.to_csv(args.outfile, sep="\t", index=False)


if __name__ == "__main__":
    main()
