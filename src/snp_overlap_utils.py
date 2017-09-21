import glob
import os
import sys

import pandas as pd

def calculate_rsquared(command, snp_list, imputation_path, r2_cutoff):
    """
    This function is designed to run the plink command for each GWAS snp

    :param command: the string command to executed via bash console
    :param snp_list: The dataframe of GWAS snps
    :param imputation_path: the path to locate the imputed snp data
    :param r2_cutoff: the cutoff r_squared value for plink to use
    :return: None but plink is executed
    """
    for index, row in snp_list.iterrows():
        return_code = os.system(command.format(imputation_path, row["Chr"], row["SNP"], r2_cutoff, row["SNP"]))
        if return_code:
            sys.stderr.write("SNP {} caused an error.\n".format(row["SNP"]))
    return


def calculate_overlap(command, snp_file, peak_df, overlap_file, window_size):
    """
    This function is designed to run the bedtools command that will determine
    what snps overlap which peak.

    :param command: the string command to executed via bash console
    :param snp_file: the data frame of imputed snps
    :param peak_df: the peak description dataframe that contains the peak file locations
    :param overlap_file: the name of the overlap file for bedtools to write to
    :param window_size: the size of the window buffer
    :return: None but bedtools is located
    """
    for index, row in peak_df.iterrows():
        return_code = os.system(command.format(row["File_Path"], window_size, snp_file, row["File_Path"], row["Tissue"], overlap_file))
        if return_code:
            sys.stderr.write("Peak {} caused an error.\n".format(row["File_Path"]))
    return

def generate_locus_dataframe_generator(overlap_df, snp_df, peak_df, specific_locus=None, keep_rsq=True, tissue=None):
    """
    This generator function is used 9in the visualization code.
    It generates an on the fly list of membership dataframes, which is a sparse matrix
    that shows what peaks a snp overlaps

    :param overlap_df: the bedtools output in dataframe form
    :param snp_df: the imputed snp dataframe
    :param peak_df: the peak desc dataframe
    :param specific_locus: if the user wants to focus on a particular locus input a string
    :param keep_rsq: if the user wants the snps to have an LD score pass true
    :param tissue: if the user wants to focus on a particular tissue input a string
    :return: a generator object that contains a list of these dataframes
    """

    for specific_tissue in peak_df["Tissue"].unique():
        tissue_peak_df = peak_df[peak_df["Tissue"] == specific_tissue]
        tissue_overlap_df = overlap_df[overlap_df["Tissue"] == specific_tissue]

        # Skip other tissues if the user wants a specific tissue
        if tissue and specific_tissue != tissue:
            continue

        # For each locus in the snp list
        for locus in snp_df["Locus"].unique():

            # Skip loci if the user wants a specific locus
            if specific_locus and specific_locus != locus:
                continue

            rsquared_list = snp_df[snp_df["Locus"] == locus]["RSquared"]
            snp_list = snp_df[snp_df["Locus"] == locus]["SNP"]
            snp_list = list(map(lambda x: x[0] + " ({})".format(x[1]), zip(snp_list, rsquared_list))) if keep_rsq else snp_list
            membership_df = pd.DataFrame(0, index=tissue_peak_df["Label"].unique(), columns=snp_list)

            for index, row in tissue_overlap_df[tissue_overlap_df["Locus"] == locus].iterrows():
                score = tissue_peak_df[tissue_peak_df["File_Path"] == row["Peak"]]["Score"].values[0]
                snp_id = row["SNP"] + " ({})".format(row["RSquared"]) if keep_rsq else row["SNP"]
                mark_label = tissue_peak_df[tissue_peak_df["File_Path"] == row["Peak"]]["Label"]
                membership_df.loc[mark_label, snp_id] += score

            yield (specific_tissue, locus, membership_df)


def write_merged_snps(snp_df, output_file):
    """
    This function writes all of plink's r_squared output into a single file.

    :param snp_df: The GWAS snps in a dataframe
    :param output_file: the final output file name
    :return: none but a unified file is generated
    """

    final_df = pd.DataFrame([], columns=["Chr", "Start", "End", "SNP", "Locus", "Sentinel", "R2"])

    # for each snp
    for snp_file in glob.glob("plink_temp/*.ld"):
        proxy_df = pd.read_csv(snp_file, delimiter=r'\s+')
        sentinel = proxy_df["SNP_A"].unique()[0]
        
        # drop the rows with no variants named
        proxy_df = proxy_df[proxy_df["SNP_B"] != '.']
        sentinel_row = snp_df[snp_df["SNP"] == sentinel]

        proxy_df["Locus"] = sentinel_row["Locus"].values[0]
        proxy_df["Sentinel"] = sentinel
        proxy_df["Chr"] = proxy_df["CHR_B"]

        proxy_df["Start"] = proxy_df["BP_B"]
        proxy_df["End"] = proxy_df["Start"] + 1
        proxy_df["SNP"] = proxy_df["SNP_B"]

        final_df = final_df.append(proxy_df[["Chr", "Start", "End", "SNP", "Locus", "Sentinel", "R2"]])

    # Sort the dataframe for easy of use with bedtools
    final_df = final_df.sort_values(["Chr", "Start"])

    # Convert the columns into the correct datatype
    # Pandas defaults to float64
    final_df["Chr"] = final_df["Chr"].astype(int)
    final_df["Start"] = final_df["Start"].astype(int)
    final_df["End"] = final_df["End"].astype(int)

    # Add chr to the beginning of each chromosome
    final_df["Chr"] = final_df["Chr"].map(lambda x: "chr" + str(x))
    final_df.to_csv(output_file, sep='\t', header=False, index=False)

    return

def write_interesting_snps(locus_generator, peak_df, snp_df, output_file):
    """
    This function is designed to write out a summary list of relavant snps after the
    bedtools system command has finished.

    :param locus_generator: the generator function for each membership dataframe
    :param peak_df: the peak desc in dataframe form
    :param snp_df: the list of imputed snps in dataframe form
    :param output_file: the name/path of the output file
    :return: none but a file is written out
    """
    critical_columns = list(peak_df[(peak_df["Priority"] == 1)]["Label"].unique())
    final_list = pd.DataFrame([], columns=["Chr", "Start", "End", "SNP", "Locus", "RSquared", "Sentinel", "Tissue"] +  critical_columns + ["Combined"])

    for tissue, locus, locus_df in locus_generator:
        sorted_df = locus_df.T.reset_index().sort_values(critical_columns, ascending=False)
        sorted_df = sorted_df[(sorted_df[critical_columns] != 0).any(axis=1)]
        final_df = pd.merge(snp_df, sorted_df, left_on="SNP", right_on="SNP")
        final_df["Tissue"] = tissue
        final_df["Combined"] = final_df[critical_columns].sum(axis=1)
        final_list = final_list.append(final_df[["Chr", "Start", "End", "SNP", "Locus", "RSquared", "Sentinel", "Tissue"] +  critical_columns + ["Combined"]])

    # Sort the data in descending order
    final_list = final_list.sort_values(critical_columns, ascending=False)

    # Convert the types to reasonable format
    final_list["Start"] = final_list["Start"].astype(pd.np.int64)
    final_list["End"] = final_list["End"].astype(pd.np.int64)
    final_list[critical_columns] = final_list[critical_columns].astype(pd.np.int64)
    final_list["Combined"] = final_list["Combined"].astype(pd.np.int64)

    # Write out the final dataset
    final_list.to_csv(output_file, sep="\t", index=False)

    return
