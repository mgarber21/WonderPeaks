import os
import subprocess
import glob
import pandas as pd


def featurecounts(directory, User_inputs_dict=None, s=1, t="CDS", T=20,
                  stardir="starout", featurecounts_subDir="featurecounts", rerun=False):
    """
    Run featureCounts to quantify gene expression from BAM files.

    Parameters:
    - s (int): Strand specificity (0 = unstranded, 1 = stranded, 2 = reverse-stranded).
    - t (str): Feature type to count (e.g., "CDS", "exon").
    - T (int): Number of threads to use.
    - User_inputs_dict (dict): Dictionary containing user input paths and parameters.
    - stardir (str): Directory containing STAR output BAM files (default: "starout").
    - featurecounts_subDir (str): Subdirectory for featureCounts output (default: "featurecounts").

    Returns:
    - int: Always returns 0 after successfully running the command.
    """
    if User_inputs_dict is None:
        raise ValueError("User_inputs_dict cannot be None")

    # Create the output directory for featureCounts results
    featurecountsDir = os.path.join(directory, featurecounts_subDir)
    os.makedirs(featurecountsDir, exist_ok=True)

    # Get the genome annotation file (GFF) and grouping key (g)
    GFF = User_inputs_dict["genome_annotations_path"]
    g = User_inputs_dict["g"]

    # Handle transcript length filters
    include_transcripts_above = User_inputs_dict.get("include_transcripts_above")
    include_transcripts_below = User_inputs_dict.get("include_transcripts_below")

    # Set input and output file paths
    filter_flag = ""
    input_files = os.path.join(directory, stardir, "*out.bam")
    output_file = os.path.join(featurecountsDir, f"{t}_STAR_counts.txt")

    if include_transcripts_above:
        filter_flag = f"above{include_transcripts_above}"
    elif include_transcripts_below:
        filter_flag = f"below{include_transcripts_below}"

    if filter_flag:
        input_files = os.path.join(directory, stardir, f"*{filter_flag}.bam")
        output_file = os.path.join(featurecountsDir, f"{t}_{filter_flag}_STAR_counts.txt")
    input_files = glob.glob(input_files)


    if os.path.isfile(output_file) and not rerun:
        return output_file
    else:
        # Run the featureCounts command
        featurecounts_cmd(T, t, g, s, GFF, output_file, input_files)
        return output_file


def featurecounts_cmd(T, t, g, s, GFF, output_file, input_files):
    """
    Construct and execute the featureCounts command.

    Parameters:
    - T (int): Number of threads to use.
    - t (str): Feature type to count (e.g., "gene", "exon").
    - g (str): Grouping key for featureCounts (e.g., gene ID).
    - s (int): Strand specificity (0 = unstranded, 1 = stranded, 2 = reverse-stranded).
    - GFF (str): Path to the genome annotation file.
    - output_file (str): Path to the output file.
    - input_files (str): Pattern for input BAM files.

    Returns:
    None
    """
    # Construct the featureCounts command
    command_list= [
        "-T", str(T),
        "-t", t,
        "--extraAttributes", "gene_biotype",
        "-g", g,
        "-s", str(s),
        "-a", GFF,
        "-o", output_file
    ]

    command = ["featureCounts"] + command_list + input_files
    # return command
    print("Running command:", " ".join(command))
    subprocess.run(command, check=True)

def featurecounts_biotype(featurecounts_file, gene_biotypes, t, PeakStream =False):
    
    if isinstance(gene_biotypes, str):
        gene_biotypes= [i.strip(" ") for i in gene_biotypes.split(";")]
    
    featurecounts_  = pd.read_csv(featurecounts_file, 
                            sep ="\t", comment = "#")

    featurecounts_gene_biotypes = featurecounts_[featurecounts_["gene_biotype"].isin(gene_biotypes)].copy()
    featurecounts_gene_biotypes.sort_values(["Chr", "Start"], inplace=True)
    gene_biotypes_label = "-".join(gene_biotypes)
    featurecounts_file_no_ext = os.path.splitext(featurecounts_file)[0]
    featurecounts_file_biotype = f"{featurecounts_file_no_ext}-{gene_biotypes_label}.txt"
    featurecounts_gene_biotypes.to_csv(featurecounts_file_biotype, index=False)
    
    return featurecounts_file_biotype
