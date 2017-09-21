import argparse

from snp_overlap_utils import *

# Set up the command line arguments
parser = argparse.ArgumentParser(description='This program is designed to find SNPS that overlap called NGS peaks.')
subparsers = parser.add_subparsers(help='sub-command help', dest='command')

# Set up Plink command
plink_parser = subparsers.add_parser('Plink', help='Use this to run Plink.')
plink_parser.add_argument('--imputation_path', '-i', type=str, help='Path that points to the imputed snps.', required=True)
plink_parser.add_argument('--r_squared', '-r', type=float, help='[Optional] Specifies an r_squared cutoff value.', default=0.4)
plink_parser.add_argument('--snp_file', '-s',  type=str, help='Path to the file containing the list of snps. Can already contain proxy snps or use'
                                                        ' plink to calculate r2 to find the proxy snps.', required=True)
plink_parser.add_argument('--output_file', '-o', type=str, help='[Optional] Path to the merged snp output file.', default="overlap_files/imputed_snps.bed")

# Set up Bedtools command
bedtools_parser = subparsers.add_parser('Bedtools', help='Use this to Find SNPS that overlap peaks.')
bedtools_parser.add_argument('--snp_file', '-s',  type=str, help='[Optional] Path to the file containing the list of snps. Can already contain proxy snps or use'
                                                        ' plink to calculate r2 to find the proxy snps.', default="overlap_files/imputed_snps.bed")
bedtools_parser.add_argument("--peak_file", '-p', type=str, help='Path to file that contains a directory of all peaks used by this program.', required=True)
bedtools_parser.add_argument('--overlap_output', '-o', type=str, help='[Optional] Path Name of the overlap file for bedtools output.', default="overlap_files/snp_overlap_new.bed")
bedtools_parser.add_argument('--window_size', '-w', type=int, help='[Optional] Specify the bp buffer for each called peak', default=50)
bedtools_parser.add_argument('--snp_output', '-so', type=str, help="[Optional] Path to write out list of overlapping snps", default="overlap_files/final_results.tsv")
bedtools_parser.add_argument('--locus', '-l', type=str, help='[Optional] Specify a single locus to write for snp output.')
bedtools_parser.add_argument('--tissue', '-t', type=str, help='[Optional] Specify what tissue to use for writing snp output.')

# Set up Heatmap command
heatmap_parser = subparsers.add_parser('Heatmap', help='Use this to generate figures.')
heatmap_parser.add_argument('--snp_file', '-s',  type=str, help='[Optional] Path to the file containing the list of snps. Can already contain proxy snps or use'
                                                        ' plink to calculate r2 to find the proxy snps.', default="overlap_files/imputed_snps.bed")
heatmap_parser.add_argument("--peak_file", '-p', type=str, help='Path to file that contains a directory of all peaks used by this program.', required=True)
heatmap_parser.add_argument('--overlap_file', '-o', type=str, help='[Optional] Name the overlap file for bedtool output.', default="overlap_files/snp_overlap_new.bed")
heatmap_parser.add_argument('--heatmap_folder', '-m', type=str, help='[Optional] Specifies the location to store the generated heatmap files.', default="overlap_files/heatmap")
heatmap_parser.add_argument('--show_plot', '-sh', help='[Optional] Have the script show an interactive heatmap graph.', action="store_true")
heatmap_parser.add_argument('--locus', '-l', type=str, help='[Optional] Specify a single locus to graph.')
heatmap_parser.add_argument('--color', '-c', type=str, help='[Optional] Path to specify the color scheme.', default='color_scheme.csv')
heatmap_parser.add_argument('--tissue', '-t', type=str, help='[Optional] Specify what tissue to use.', default='Fat')
heatmap_parser.add_argument('--cube', '-cu', help='[Optional] Specify if the program should print 3D.', action="store_true")

# Parse the arguments
args = parser.parse_args()

# Set up the directory for this program
if not os.path.exists("overlap_files"):
    os.mkdir("overlap_files")
    os.chdir("overlap_files")
    os.mkdir("heatmap")
    os.chdir("..")

# If the Plink command is used
if args.command == 'Plink':

    # Output for the user to monitor progress
    print("Gathering Proxy Snps")

    # Set up the command and read the snp files
    plink_command = "plink --vcf {}/{}.vcf.gz --ld-snp {} --ld-window-r2 {} --out {}_LD_block --double-id --r2 inter-chr with-freqs --memory 6000"
    snp_df = pd.read_csv(args.snp_file, delimiter=r'\s+',engine='python')

    # Create a folder to write the plink files out
    # Plus easy deletion once done
    if not os.path.exists("plink_temp"):
        os.mkdir("plink_temp")

    # Contain the massive files that are going to be generated
    os.chdir("plink_temp")

    # Run the plink command
    calculate_rsquared(plink_command, snp_df, args.imputation_path, args.r_squared)

    #Change back to parent directory
    os.chdir("..")

    # Merge the individual ld files
    write_merged_snps(snp_df, args.output_file)

# If calculating overlap
if args.command == 'Bedtools':

    # Output for the user to monitor progress
    print("Calculating SNP Overlap")

    # Set up the system command
    peak_df = pd.read_csv(args.peak_file, delimiter=r'[\t|,]', engine="python")
    bedtools_command = "cut -d$'\\t' -f 1,2,3 {}"
    bedtools_command = bedtools_command + " | bedtools window -w {} -a stdin -b {}"
    bedtools_command =  bedtools_command + " | awk '{{print $0 \"\t{}\t{}\"}}' >> {}"

    # Run the command
    calculate_overlap(bedtools_command, args.snp_file, peak_df, args.overlap_output, args.window_size)

    # Write out the snps that overlap peaks
    print("Writing Relevant SNPS")

    # Name of the columns pandas will use
    snp_names = ["Chr", "Start", "End", "SNP", "Locus", "Sentinel", "RSquared"]
    overlap_names = ["ChrB", "ChrB_Start", "ChrB_End", "ChrA", "ChrA_Start", "ChrA_End", "SNP", "Locus", "Sentinel",
                     "RSquared", "Peak", "Tissue"]

    # Open the files
    snp_df = pd.read_csv(args.snp_file, names=snp_names, sep="\t")
    overlap_df = pd.read_csv(args.overlap_output, names=overlap_names, sep="\t", index_col=False, header=None).drop_duplicates()
    overlap_df["ChrB_Start"] = overlap_df["ChrB_Start"].apply(lambda x: x - args.window_size)
    overlap_df["ChrB_End"] = overlap_df["ChrB_End"].apply(lambda x: x + args.window_size)
    overlap_df.to_csv(args.overlap_output, header=False, index=False, sep="\t")

    # Create generator and write the snps out
    locus_generator = generate_locus_dataframe_generator(overlap_df, snp_df, peak_df, specific_locus=args.locus, keep_rsq=False, tissue=args.tissue)
    write_interesting_snps(locus_generator, peak_df, snp_df, args.snp_output)

# If the heatmap command is used
if args.command == 'Heatmap':

    # Conditional import for heatmap command
    from colour import Color

    # Output for the user to monitor progress
    print("Generating HeatMap")

    # Get the peak definitions file
    peak_df = pd.read_csv(args.peak_file, delimiter=r'[\t|,]', engine="python")

    # Name of the columns pandas will use
    snp_names = ["Chr", "Start", "End", "SNP", "Locus", "Sentinel", "RSquared"]
    overlap_names = ["ChrB", "ChrB_Start", "ChrB_End","ChrA", "ChrA_Start", "ChrA_End", "SNP", "Locus", "Sentinel", "RSquared",
                     "Peak", "Tissue"]

    # Open the files
    snp_df = pd.read_csv(args.snp_file, names=snp_names, sep="\t")
    overlap_df = pd.read_csv(args.overlap_file, names=overlap_names, sep="\t", index_col=False, header=None)

    # Sort the Color map and cast each string into a Color object
    flatui = map(lambda x: Color(x), pd.read_csv(args.color, delimiter=r'[\t|,]', engine="python").sort_values("Rank")["Colors"])

    # Create the generator and generate a heatmap for each locus
    if args.cube:
        # conditional import for 3D
        # This is a literal use for saying the 3D vis is not implemented
        raise NotImplementedError()

        #from snp_overlap_vis_3D import *
        #generate_heatmap_cube(generate_locus_dataframe_generator, overlap_df, snp_df, peak_df, args.heatmap_folder, flatui, specific_locus=args.locus, show_plots=args.show_plot)
    else:
        # conditional import for 2D
        from snp_overlap_vis import *

        locus_generator = generate_locus_dataframe_generator(overlap_df, snp_df, peak_df, specific_locus=args.locus, keep_rsq=True, tissue=args.tissue)
        generate_heatmap(locus_generator, snp_df, peak_df, args.heatmap_folder, flatui, show_plots=args.show_plot)
