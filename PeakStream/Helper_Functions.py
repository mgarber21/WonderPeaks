import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import hypergeom
import csv
from scipy.misc import electrocardiogram
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import seaborn as sns
import math
from Bio.SeqIO.FastaIO import SimpleFastaParser


def PeakStream_Direcory(directory):
    PeakStream_directory = os.path.join(directory, "PeakStream")
    os.makedirs(PeakStream_directory, exist_ok=True)
    return PeakStream_directory


def q25(x):
    # 25th Percentile
    return x.quantile(0.25)
def q50(x):
    # 50th Percentile
    return x.quantile(0.5)
def q75(x):
    # 75th Percentile
    return x.quantile(0.75)
def q90(x):
    # 90th Percentile
    return x.quantile(0.9)


def stream(row, gene_buffer  = 600, step = 100, maxUTR = 5000):
    """
    Link peaks to genes within defined boundaries using strand-specific logic.

    Parameters:
    - row (pd.Series): Row representing a genomic region and peaks.
    - gene_buffer (int): Distance to exclude from gene boundaries.
    - step (int): Step size to extend search boundaries.
    - maxUTR (int): Maximum allowable UTR length.

    Returns:
    - dict: Links of peaks to genomic features.
    """

    strand_multiplier = dict(zip(["-", "+"], [-1, 1]))
    strand_start = dict(zip(["-", "+"], ["end", "start"]))
    strand_end = dict(zip(["-", "+"], ["start", "end"]))
    strand = row["strand"]
    strand_next_gene_end = dict(zip(["-", "+"], ["start_preceding_gene", "end_next_gene"]))
    strand_next_gene_length = dict(zip(["-", "+"], ["preceding_gene_length", "next_gene_length"]))

    iter_max = int(maxUTR/step)

    o_end = row[strand_end[strand]] - strand_multiplier[strand]*min(gene_buffer, row["gene_length"])
    o_end_next_gene = row[strand_next_gene_end[strand]] - strand_multiplier[strand]*min(gene_buffer, row[strand_next_gene_length[strand]])

    dict_peak_link = dict()
    for i in range(iter_max):
        n_end = o_end + strand_multiplier[strand]*step*i
        peak = row['peak_location']
        
        # if peak fall within the bound of the next gene up (+) or downstream (-)
        if strand_multiplier[strand]*peak > strand_multiplier[strand]*o_end_next_gene:
            break
        
        if strand_multiplier[strand]*peak > strand_multiplier[strand]*o_end and strand_multiplier[strand]*peak < strand_multiplier[strand]*n_end:
            dict_peak_link[peak] = [row["seqname"], strand, row["start"], row["end"], row["score_peak"]]

    return dict_peak_link

def load_coordinates(coordinate_file):
    """
    Load and preprocess gene coordinates from a GTF file.

    Returns:
    - pd.DataFrame: Dataframe of gene coordinates with additional metadata.
    """

    gtf_data = pd.read_csv(coordinate_file, index_col=False)
    
    # Attempt to filter by "gene"
    if "gene" in gtf_data["feature"].unique():
        gtf_data = gtf_data[gtf_data["feature"] == "gene"]
    elif "exon" in gtf_data["feature"].unique():
        gtf_data = gtf_data[gtf_data["feature"] == "exon"]
    else:
        raise ValueError("Neither 'gene' nor 'exon' features found in the annotations file.")
    
    gtf_data.sort_values(["seqname", "start","end"], inplace = True)

    gtf_data["start_next_gene"] = gtf_data.groupby(["seqname", "strand"])["start"].shift(periods =-1)
    gtf_data["strand_next_gene"] = gtf_data.groupby(["seqname", "strand"])["strand"].shift(periods =-1)
    gtf_data["end_next_gene"] = gtf_data.groupby(["seqname", "strand"])["end"].shift(periods =-1)


    gtf_data["end_preceding_gene"] = gtf_data.groupby(["seqname", "strand"])["end"].shift(periods =1)
    gtf_data["strand_preceding_gene"] = gtf_data.groupby(["seqname", "strand"])["strand"].shift(periods =1)
    gtf_data["start_preceding_gene"] = gtf_data.groupby(["seqname", "strand"])["start"].shift(periods =1)


    # Apply the function and unpack the results into two new columns
    gtf_data[["TandemGene", "dist_to_TandemGene"]] = gtf_data.apply(
        lambda row: pd.Series(Tandem_gene(row)), axis=1
    )
    return gtf_data

def rename_attributes(att, named_attribute = None):

    if "=" in att:
        # Assume GFF formatting
        dict_attributes = {
            i.split("=")[0]: i.split("=")[1]
            for i in att.split(";")
            if len(i.split("=")) > 1
        }
        
        # Try to extract the ID, Name, or Parent in that order
        ID = None


        if named_attribute:
            try:
                ID = dict_attributes.get(named_attribute)
                return ID
            except Exception as e:
                print(f"Error extracting {named_attribute} from annotations file: {e}")
        try:
            ID = dict_attributes.get("ID") or dict_attributes.get("Name") or dict_attributes.get("Parent")
        except Exception as e:
            print(f"Error extracting ID from annotations file: {e}")


    else:
        # assume this means GTF formatting
        dict_attributes = {i.split(" ")[0].strip('''"'''):i.split(" ")[1].strip('''"''').strip(''';"''') for i in att.split("; ") if len(i.split(" ")) > 1}

        # Try to extract the gene_id, transcript_id, or gene_name in that order
        ID = None
        if named_attribute:
            try:
                ID = dict_attributes.get(named_attribute)
                return ID
            except Exception as e:
                print(f"Error extracting {named_attribute} from annotations file: {e}")
        

        try:
            ID = dict_attributes.get("gene_id") or dict_attributes.get("transcript_id") or dict_attributes.get("gene_name")
        except Exception as e:
            print(f"Error extracting ID from annotations file: {e}")

    return ID

def Tandem_gene(row):
    """
    Determines if the gene is in tandem with the adjacent gene 
    and calculates the distance to the tandem gene.

    Parameters:
    - row (pd.Series): A row of the GTF DataFrame.

    Returns:
    - tuple: (bool, float or np.nan)
        - True and distance if the gene is in tandem.
        - False and NaN otherwise.
    """
    if row["strand"] == "+" and row["strand_next_gene"] == "+":
        dist_to_next = row["start_next_gene"] - row["end"]
        return (True, dist_to_next)
    elif row["strand"] == "-" and row["strand_preceding_gene"] == "-":
        dist_to_preceding = row["start"] - row["end_preceding_gene"]
        return (True, dist_to_preceding)
    else:
        return (False, np.nan)  # Explicit return statement for non-tandem cases


def outdir(PS_directory, new_coordinates_file_directory):
        if not new_coordinates_file_directory:
            # infer the directory name from the input folder
            
            new_coordinates_file_directory = PS_directory

        os.makedirs(new_coordinates_file_directory, exist_ok=True)
        return new_coordinates_file_directory

    
def prefix(coordinate_file, coordinates_file_prefix):
    if not coordinates_file_prefix:
        coordinates_file_prefix = os.path.basename(coordinate_file).split(".")[0]
    return coordinates_file_prefix

def suffix(PS_directory, coordinates_file_suffix):
    if not coordinates_file_suffix:
        coordinates_file_suffix = os.path.basename(os.path.dirname(PS_directory))
    return coordinates_file_suffix

def get_basename(filename):
    basename = os.path.basename(os.path.splitext(filename)[0]\
        .split(".fastq")[0]\
            .split("_fastp_output")[0])
    
    return basename

def metadata_upload(directory,designfactor,
                    meta_delim_whitespace = False,
                    meta_index_col = None):
        
    meta_file_path = os.path.join(directory, "NGS_user_metadata.csv")
    metadata = pd.read_csv(meta_file_path,
                delim_whitespace=meta_delim_whitespace, 
                index_col=meta_index_col)\
                .rename_axis(mapper= None,axis=0)

    metadata["basename"] =  metadata.apply(lambda row:
                        get_basename(row["rawdata_filename"]), axis =1)


    metadata["StrainCondition"] = metadata.apply(lambda row: f"{row['Strain']}_{row['Condition']}", axis=1)
    metadata["SampleCondition"] = metadata.apply(lambda row: f"{row['Sample']}_{row['Condition']}", axis=1)
    metadata["SampleStrain"] = metadata.apply(lambda row: f"{row['Strain']}_{row['Sample']}", axis=1)

    return metadata