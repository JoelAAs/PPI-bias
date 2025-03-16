import pandas as pd

def set_relative_frequency(infile, outfile):
    df_intensities = pd.read_csv(infile, sep=",")
    df_intensities = df_intensities.set_index("Gene names")
    relative_freq = df_intensities.div(
        df_intensities.sum(
            axis=0,
            skipna=True)).mean(axis=1)
    relative_freq.to_csv(
        outfile,
        sep="\t",
        index_label=["gene_name", "relative_frequence"]
    )

set_relative_frequency(
    "intensities_wide_selected_N04547_M07444.csv",
    "relative_freq_wide_selected_N04547_M07444.csv"
    )

#set_relative_frequency(
#    "intensities_wide_selected_N07444_M04547.csv",
#    "relative_freq_wide_selected_N07444_M04547.csv"
#    )