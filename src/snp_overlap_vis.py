from matplotlib.colors import ListedColormap

import matplotlib
matplotlib.use('Agg')

import seaborn as sns

def generate_heatmap(locus_generator, snp_df, peak_df, heatmap_folder, color_scheme, show_plots=True):
    """
    This function is designed to generate an overlap heatmap.

    :param locus_generator: the generator for each membership dataframe
    :param snp_df: the imputed snps in dataframe form
    :param peak_df: the peak desc file in a dataframe
    :param heatmap_folder: the folder that will contain all the generated heatmaps
    :param color_scheme: the color scheme list for the heatmap
    :param show_plots: a boolean to interactively show each generated heatmap usually safe to keep false
    :return: no object just a lot of heatmaps are generated
    """
    critical_columns = list(peak_df[peak_df["Priority"] == 1]["Label"].unique())

    for tissue, locus, locus_df in locus_generator:
        display_df = locus_df

        display_df = display_df.T.sort_values(critical_columns, ascending=False)
        display_df = display_df[(display_df[critical_columns] != 0).any(axis=1)]

        # Skip df if empty
        if display_df.empty:
            print "Locus: {} in Tissue: {} doesn't have snps. Skipping!!!".format(locus, tissue)
            continue

        # Change the figure size
        sns.set_context("paper")
        sns.plt.figure(figsize=(10, 6))

        ax = sns.heatmap(display_df,
                         cmap=ListedColormap(map(lambda x: x.rgb, color_scheme)), linewidths=0.5, cbar=True, square=False,
                         vmin=-1)

        ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=10, fontsize='x-small')
        ax.set_yticklabels(ax.yaxis.get_majorticklabels(), rotation=0, fontsize=7)
        ax.set_title("Peak Overlap Heatmap ({})\n Total # SNPS: {} \n Sentinel SNP: {}".format(
            locus, snp_df[snp_df["Locus"] == locus].shape[0],
            snp_df[snp_df["Locus"] == locus]["Sentinel"].values[0])
        )

        sns.plt.savefig("{}/{}.png".format(heatmap_folder, locus))

        if show_plots:
            sns.plt.show()

        sns.plt.close()
