"""
Implements peak detection on bedgraph data using the first derivative and integrates gene annotations.

"""

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

from PeakStream.Helper_Functions import *

class WonderPeaks_4UTRs():
    """
    Process bedgraph files and detect transcript peaks using the first derivative.
    
    """

    def __init__(self, directory, coordinate_file, n = 10):

        """
        Initialize and load genome sequences and bedgraph files for peak calling.

        Parameters:
        - directory (str): path to directory with data
        - coordinate_file (str): path to genome coordinate file
        - n (int): steps used to calculate peaks

        """

        
        self.gtf_data  = load_coordinates(coordinate_file = coordinate_file)
        self.directory = directory
        self.PS_directory = PeakStream_Direcory(self.directory)
        self.n = n
        self.bed_directory = os.path.join(self.directory, "bedgraphout")

        

    def load_files(self):
        """
        Initialize bedgraph files for peak calling.

        returns: 
        - dict: dictionary of files in your study
        
        """

        strands = ["fwd", "rev"]
        files = []
        bedfiles_dict  = {}
        for strand in strands:
            
            files = []
            for file in [file for file in os.listdir(self.bed_directory) if strand in file]:
                files.append( os.path.join(self.bed_directory,file))
                files = sorted(files)
            bedfiles_dict[strand] = files

        return bedfiles_dict
    
    def load_data(self, bedfiles_dict):

        """
        Load bedgraph files as dataframes for peak calling.

        returns: 
        - pd.DataFrame: DataFrame of your concatenated begraph files
        
        """

        # bedfiles_dict = self.load_files()

        # global is set for all functions in this page
        global df

        df = pd.DataFrame()
        for strand, files in bedfiles_dict.items():
            for file in files:
                df_bed = pd.read_csv(file, sep = "\t", header = None)
                df_bed.rename(columns=dict(zip(df_bed.columns, ["chr", "start", "stop", "score"])), inplace =True)


                df_bed["file"] = file
                df_bed["dx"], df_bed["dy"] = df_bed.groupby(['file', 'chr'])['start'].diff().fillna(20), df_bed.groupby(['file', 'chr'])['score'].diff().fillna(0)

                df_bed["1st_derivative"] = df_bed["dy"]/df_bed["dx"]
                df_bed["1st_derivative"] = df_bed["1st_derivative"].interpolate(method = 'linear')

                df_bed.fillna(0,inplace = True)
                df_bed.dropna(how = "any", inplace=True)
                df_bed["strand"] = strand
                df = pd.concat([df,df_bed])
            

        df.replace([np.inf, -np.inf], 0, inplace=True)


        self.df_summary = df.groupby("file")["score"].agg(["median","mean", q25, q50, q75, q90]).reset_index()
        self.df_summary.to_csv(os.path.join(self.PS_directory, "bedgraph_summary.csv"), index = False)


        return df
    
    def first_derivative(self, df):
        """
        Take the first derivative of each file to detect peaks (location of transcripts)

        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files
        - pd.DataFrame: DataFrame of your concatenated begraph files
        
        """

        # n =10  # number of points to be checked before and after the peak. Using 10 bc its half the binsize (20).    
        
        # df = self.load_data()

        df_peaks = pd.DataFrame()

        for chr in df["chr"].unique():
            df_chr = df[df["chr"] == chr].reset_index()

            # Find local peaks

            df_chr['min'] = df_chr.iloc[argrelextrema(df_chr["1st_derivative"].values, np.less_equal,
                                order=self.n)[0]]['1st_derivative']
            df_chr['max'] = df_chr.iloc[argrelextrema(df_chr["1st_derivative"].values, np.greater_equal,
                                order=self.n)[0]]['1st_derivative']

            
            df_peaks_chr = df_chr[
                ((df_chr["max"] > 0.005) | (df_chr["min"] < -0.005))
                    ]

            df_peaks = pd.concat([df_peaks, df_peaks_chr]).reset_index(drop=True)



        df_peaks = df_peaks[list(df.columns) + [ "max", "min"]]



        return df_peaks, df    
    

    def median_score_cutoff(self, df_peaks_compared_merged):

        """
        Uses the median score of the bedgraph file (from the summary file) to determine whether a peak is True
        by setting this value as the threshold.

        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files with a score above a certain threshold
        
        """

        df_peaks_compared_merged_thresh = pd.merge(df_peaks_compared_merged, self.df_summary, on = "file")

        df_peaks_compared_merged_thresh= df_peaks_compared_merged_thresh[
            df_peaks_compared_merged_thresh["score_peak"] > df_peaks_compared_merged_thresh["median"]
            ]

        return  df_peaks_compared_merged_thresh

    def clean_call_peaks(self, peaks):
        """
        Applies functions to rename columns and cleanup data frames for next steps
        Applies a threshold (500bp) for the width of a True peak

        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files above a certain width threshold
        
        """
        peaks['nucleotide'] = peaks.apply(lambda row: "A" if row['strand'] == "fwd" else "T", axis= 1)
        peaks['strand_direction'] = peaks['strand']
        peaks['strand'] = peaks.apply(lambda row: "+" if row['strand'] == "fwd" else "-", axis= 1)
        peaks.rename(columns={"chr":"seqname"}, inplace =True)
        
        # if the length of the peak start and end is > 500, it is not a true peak
        peaks[peaks["peak_diff"] < 500]

        return peaks
    


    def normalize_peak_location(self, peaks):

        """
        Applies functions to rename columns and cleanup data frames for next steps
        Applies a threshold (500bp) for the width of a True peak

        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files above a certain width threshold
        
        """
        peaks.sort_values(["seqname", "peak_location", "score"], inplace=True)
        
        peaks["peak_location_actual"] = peaks["peak_location"]
        
      #   try another method to get peak locations
              # set peak locations within 30bp of eachother to eachother

        # Set tolerance
        print("Centering peaks...")
        
        # set peak locations within 30bp of eachother to eachother

        peaks.sort_values(["seqname", "peak_location"], inplace=True)
        peaks["peak_location_ranges_diff_by_rep"] =peaks.groupby("seqname")["peak_location"].diff()
        peaks['peak_location_ranges_diff_by_rep'].fillna(0, inplace=True)

        peaks['peak_location_ranges_diff_by_rep'].replace(to_replace=0, method='ffill', inplace=True)

        peaks["peak_location_LeftAlign"] = peaks.apply(lambda row : row["peak_location"]-row["peak_location_ranges_diff_by_rep"] \
                                                                     if abs(row["peak_location_ranges_diff_by_rep"])<= 100 \
                                                                        else row["peak_location"], 
                                                                        axis =1)
        peaks["peak_location_Centered"] = peaks.groupby(["seqname","strand",  "peak_location_LeftAlign"])["peak_location"].transform('mean')

        peaks["peak_location"]=peaks["peak_location_Centered"] 

        return peaks

    def count_peak_occurances(self, peaks):

        """
        Applies function to count the number of peaks across bedgraph files
        and sets a threshold to only include peaks that occur in more than 2 files


        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files above a certain count threshold,
        
        """

        # count the number of occurances of a peak_location
        # peak must occur more than twice to be considered for linking

        # peak_location_counts = pd.DataFrame(peaks[["seqname","strand", "peak_location"]].value_counts()).rename(columns = {0:"peak_location_counts"}).reset_index()
        peaks["peak_location_counts"] = peaks.groupby(["seqname","strand",  "peak_location"])["peak_location"].transform('count')
        # peaks = pd.merge(peaks, peak_location_counts, on =["seqname", "strand", "peak_location"])
        count_centeredpeaks = ((peaks["peak_location_counts"] > 2))
        peaks = peaks[count_centeredpeaks]

        return peaks
    
    def finalize_call_peaks(self, peaks):
        """
        Applies functions to apply thresholds to data


        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files above a thresholds,
        
        """
        peaks = self.clean_call_peaks(peaks)
        peaks = self.normalize_peak_location(peaks)
        peaks = self.count_peak_occurances(peaks)

        return peaks

    def call_peaks(self, median_score_cutoff_ = True):
        """
        Applies first derivative to bedgraph files to calculate the occurance of a peak based on passing through y=0
        by taking into account the location of neig peaks in the 1st derivative

        Parameters:
        - median_score_cutoff_ (bool): whether to use the median score of the bedgraph file as a cutoff value for peak detection
        
        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files above a certain count threshold,
        
        """

        print("Running WonderPeaks for UTR assingment...")
        pd.options.mode.chained_assignment = None
        
        # run funtions to call peaks
        bedfiles_dict = self.load_files()
        df = self.load_data(bedfiles_dict)
        df_peaks, df = self.first_derivative(df)
        

        df_peaks["next_min"] = df_peaks.groupby(["file","chr"])["min"].shift(periods =-1)
        df_peaks["previous_max"] = df_peaks.groupby(["file","chr"])["max"].shift(periods =1)
        df_peaks["previous_stop"] = df_peaks.groupby(["file","chr"])["stop"].shift(periods =1)
        df_peaks.reset_index(inplace=True, drop =True)
        df_peaks_compared = df_peaks[
            ((~df_peaks["next_min"].isna()) & (~df_peaks["max"].isna()))
            |
            ((~df_peaks["previous_max"].isna()) & (~df_peaks["min"].isna()))
            |
            ((~df_peaks["previous_max"].isna()) & (~df_peaks["min"].isna()) & (~df_peaks["next_min"].isna()) & (~df_peaks["max"].isna()))
            ]



        df_peaks_compared.reset_index(inplace=True, drop =True)

        df_peaks_compared["peak_diff"] = df_peaks_compared.groupby(["file","chr"])["start"].diff()
        df_peaks_compared.dropna(subset = ["min"], inplace=True)
        df_peaks_compared["peak_location"] = df_peaks_compared["stop"] - df_peaks_compared["peak_diff"] + 2*self.n
        df_peaks_compared_merged = pd.merge(df_peaks_compared, df, left_on =["chr", "file", "peak_location"], 
                                                right_on = ["chr","file","start"], how = "inner", suffixes  =("", "_peak"))

        if median_score_cutoff_:
            df_peaks_compared_merged_thresh = self.median_score_cutoff(df_peaks_compared_merged)
            return self.finalize_call_peaks(peaks = df_peaks_compared_merged_thresh)
        

        return self.finalize_call_peaks(peaks = df_peaks_compared_merged)