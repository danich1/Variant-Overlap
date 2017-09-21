import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def generate_heatmap_cube(locus_generator, overlap_df, snp_df, peak_df, heatmap_folder, color_scheme, specific_locus=None, show_plots=True):
    """
    This function is under maintence as in I ran out of time to finish this. Nevertheless.
    the goal of this function is to show a 3D version of the overlap memberships. Ideally,
    this will be an active feature once the data for more tissues becomes available.

    :param locus_generator: the generator for each membership dataframe
    :param overlap_df: the bedtools output in dataframe form
    :param snp_df: the imputed snps in dataframe form
    :param peak_df: the peak desc file in a dataframe
    :param heatmap_folder: the folder that will contain all the generated heatmaps
    :param color_scheme: the color scheme list for the heatmap
    :param specific_locus: if the user wants to focus on a particular locus input a string
    :param show_plots: a boolean to interactively show each generated heatmap usually safe to keep false
    :return: none but a 3D image will be generated
    """

    values = map(lambda x: list(x.rgb) + [1], color_scheme)
    values[3][3] = 0.2
    tissues = peak_df["Tissue"].unique()

    if specific_locus:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        for i, data in enumerate(locus_generator(overlap_df, snp_df, peak_df, specific_locus=specific_locus, keep_rsq=False)):
            tissue, current_locus, locus_df = data
            display_df = locus_df
            if i == 0:
                critical_columns = list(peak_df[(peak_df["Tissue"] == tissue) & (peak_df["Priority"] == 1)]["Label"].unique())
                display_df = display_df.T.sort_values(critical_columns, ascending=True)
                display_df = display_df[(display_df[critical_columns] != 0).any(axis=1)]
                display_df = display_df.T
                selected_snps = list(display_df.columns)

                mark_length = len(peak_df["Label"].unique())
                snp_length = display_df.shape[1]

                X = np.arange(0, mark_length, 1)
                Y = np.arange(0, len(tissues), 1)
                Z = np.arange(0, snp_length, 1)
                X, Y, Z = np.meshgrid(X, Y, Z)
            else:
                display_df = locus_df[selected_snps]

            surface_map = np.zeros((mark_length, snp_length, 4))

            display_df = display_df + 3
            for z, snp in enumerate(display_df.columns):
                for x, mark in enumerate(display_df.index):
                    location = display_df.loc[mark][snp]
                    surface_map[x][z] = values[location]
            surf = ax.plot_surface(X[i], Y[i], Z[i], facecolors=surface_map)

        ax.set_xticklabels(list(display_df.index))
        ax.set_yticks(np.arange(0, 8, 1))
        ax.set_yticklabels(tissues)
        ax.set_zticks(np.arange(0, display_df.shape[1], 1))
        ax.set_zticklabels(list(display_df.columns))

        plt.savefig("{}/{}_3D.png".format(heatmap_folder, locus))

    else:
        for locus in snp_df["Locus"].unique():
            fig = plt.figure()
            fig.subplots_adjust(bottom=0.2)
            ax = fig.gca(projection='3d')

            # TODO create a function that will return a list of dataframes and the important snps
            for i, data in enumerate(locus_generator(overlap_df, snp_df, peak_df, specific_locus=locus, keep_rsq=False)):
                tissue, current_locus, locus_df = data
                display_df = locus_df

                if i == 0:
                    critical_columns = list(peak_df[(peak_df["Tissue"] == tissue) & (peak_df["Priority"] == 1)]["Label"].unique())
                    display_df = display_df.T.sort_values(critical_columns, ascending=True)
                    display_df = display_df[(display_df[critical_columns] != 0).any(axis=1)]
                    display_df = display_df.T
                    selected_snps = list(display_df.columns)

                    mark_length = len(peak_df["Label"].unique())
                    snp_length = display_df.shape[1]

                    X = np.arange(0, mark_length, 1)
                    Y = np.arange(0, len(tissues), 1)
                    Z = np.arange(0, snp_length, 1)
                    X, Y, Z = np.meshgrid(X, Y, Z)

                else:
                    display_df = locus_df[selected_snps]

                if display_df.empty:
                    break

                surface_map = np.zeros((mark_length, snp_length, 4))
                print surface_map.shape

                display_df = display_df + 3
                for z, snp in enumerate(display_df.columns):
                    for x, mark in enumerate(display_df.index):
                        location = display_df.loc[mark][snp]
                        surface_map[x][z] = values[location]

                print Y.shape

                surf = ax.plot_surface(X[i], Y[i], Z[i], facecolors=surface_map)

            print display_df.shape
            if display_df.empty or display_df.shape[0] < 2:
                continue

            ax.set_xticklabels(list(display_df.index))
            ax.set_yticks(np.arange(0, 8, 1))
            ax.set_yticklabels(tissues)
            ax.set_zticks(np.arange(0, display_df.shape[1], 1))
            ax.set_zticklabels(list(display_df.columns))
            ax.set_title("Locus: {}".format(locus))

            plt.savefig("{}/{}_3D.png".format(heatmap_folder, locus))

            plt.clf()
            plt.close()
