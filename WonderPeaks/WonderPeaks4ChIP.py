# WonderPeaks for ChiP-seq 

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns


# from scipy.stats import hypergeom
# from scipy.misc import electrocardiogram
# from scipy.signal import find_peaks
from scipy.signal import argrelextrema

import math

from Helper_Functions import *


class WonderPeaks4ChIP():
     def __init__(self, directory, coordinate_file,
                  designfactors,slope_thresh,
                  n=10):
        """
        Initialize and load genome sequences and annotations WonderPeaks.

        Parameters:
        - coordinate_file (str): path to genome coordinate file
        - n (int): steps for 1st derivative 
        """
        
        self.directory = directory
        self.WP_directory = WonderPeaks_Direcory(self.directory)
        self.n = n
        self.bed_directory = os.path.join(self.directory, "bedgraphout", "normalizeUsingCPM")
        self.designfactors = designfactors
        self.metadata = metadata_upload(self.directory, self.designfactors)
        self.slope_thresh = slope_thresh
        self.WOnder_init = os.path.join(self.WP_directory, "WOnder_init.csv")
        self.WOnder_unfiltered_peaks = os.path.join(self.WP_directory, "WOnder_unfiltered_peaks.csv")
        
        
    
     def load_files(self):
        """
        Initialize bedgraph files for peak calling.

        returns: 
        - dict: dictionary of files in your study
        
        """
        
            
        files_list = []
        for file in os.listdir(self.bed_directory):
            files_list.append(os.path.join(self.bed_directory,file))
            files_list = sorted(files_list)
            
        
        return files_list
    
     def load_data(self, files_list):
        """
        Load bedgraph files as dataframes for peak calling.

        returns: 
        - pd.DataFrame: DataFrame of your concatenated begraph files
        
        """
        
        if os.path.isfile(self.WOnder_init):
            print("Reading in 1st derivative of trace data")
            df = pd.read_csv(self.WOnder_init)
            self.df_summary = pd.read_csv(os.path.join(self.WP_directory, "bedgraph_summary.csv"))
            return df
        
        print("Calculating 1st derivative of trace data")
        df = pd.DataFrame()
        
        for file in files_list:
           
            
            df_bed = pd.read_csv(file, sep = "\t", header = None)
            df_bed.rename(columns=dict(zip(df_bed.columns, ["seqname", "start", "stop", "score"])), inplace =True)
            df_bed_summary = df_bed.groupby(["seqname"])["score"].agg(["median","mean", q25, q50, q75, q90]).reset_index()
            
            df_bed["file"] = file
            df_bed["basename"] = get_basename(file)
            df_bed = pd.merge(df_bed,df_bed_summary, on = ["seqname"])
            
            # calculate the noramlized score by dividing the score by the chromosomes median score 
            df_bed["normalized_score"] = df_bed["score"]/df_bed["median"]

            # calculate dx and dy
            df_bed["dx"], df_bed["dy"] = df_bed.groupby(['file', 'seqname'])['start'].diff().fillna(20), df_bed.groupby(['file', 'seqname'])['score'].diff().fillna(0)
 
            # calculate 1st derivative dy/dx and interpolate
            df_bed["1st_derivative"] = df_bed["dy"]/df_bed["dx"]
            df_bed["1st_derivative"] = df_bed["1st_derivative"].interpolate(method = 'linear')

            df_bed.fillna(0,inplace = True)
            df_bed.dropna(how = "any", inplace=True)
            
            df = pd.concat([df,df_bed])
        df.replace([np.inf, -np.inf], 0, inplace=True)


        # save data for next steps and faster run times
       
        self.df_summary = df.groupby(["file", "seqname"])["score"].agg(["median","mean", q25, q50, q75, q90]).reset_index()
        self.df_summary.to_csv(os.path.join(self.WP_directory, "bedgraph_summary.csv"), index = False)
        df.to_csv(self.WOnder_init, index = False)

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

        for seqname in df["seqname"].unique():
            df_seqname = df[df["seqname"] == seqname].reset_index()

            # Find local peaks

            df_seqname['min'] = df_seqname.iloc[argrelextrema(df_seqname["1st_derivative"].values, np.less_equal,
                                order=self.n)[0]]['1st_derivative']
            df_seqname['max'] = df_seqname.iloc[argrelextrema(df_seqname["1st_derivative"].values, np.greater_equal,
                                order=self.n)[0]]['1st_derivative']

            
            df_peaks_seqname = df_seqname[
                ((df_seqname["max"] > self.slope_thresh) | (df_seqname["min"] < -self.slope_thresh))
                    ]

            df_peaks = pd.concat([df_peaks, df_peaks_seqname]).reset_index(drop=True)



        df_peaks = df_peaks[list(df.columns) + [ "max", "min"]]



        return df_peaks, df    
    

     def score_cutoff(self, df_peaks_compared_merged):

        """
        Uses the median score or percentile score of the bedgraph file (from the summary file) to determine whether a peak is True
        by setting this value as the threshold.
        

        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files with a score above a certain threshold
        
        """

        
        df_peaks_compared_merged_thresh = pd.merge(df_peaks_compared_merged, self.df_summary, on = "file")
        
        df_peaks_compared_merged_thresh= df_peaks_compared_merged_thresh[
            df_peaks_compared_merged_thresh["score_peak"] > df_peaks_compared_merged_thresh["median"]
            ]
        

        return  df_peaks_compared_merged_thresh


      
      


     def peak_shift_value(self, row, LeftAlign_bp):
      #   ensures peak will be within 120bp of next peak called
        peak_diff = max(150, row["peak_diff"])
        LeftAlign_value = LeftAlign_bp*round(peak_diff/100)
        return LeftAlign_value
     
     def shift_peak_location(self, peaks):

        """
        Applies functions to shift peak locations that are with in 80bp of one another

        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files above a certain width threshold
        
        """
        peaks.sort_values(["seqname", "peak_location", "score"], inplace=True)
        peaks["peak_location_actual"] = peaks["peak_location"]
        
      #   try another method to get peak locations
        # Set tolerance
        tolerance, broad_tolerance = 120, 300

         # Sort values for proper grouping
        peaks.sort_values(["seqname", "peak_location_actual", "score"], inplace=True)
        
      #   written with help from ChatGPT
        for tolerance, tolerance_type in zip([tolerance, broad_tolerance], ["", "_broad"]):

            # Calculate group labels with groupby and transform
            peaks[f"peak_location_group{tolerance_type}"] = peaks.groupby("seqname").apply(
                  lambda group: ((group["peak_location_actual"].diff().abs() > tolerance) #difference between the previous row
                                 ).cumsum()
               ).reset_index(level=0, drop=True)
            

            
            peaks[f"peak_location{tolerance_type}"] = peaks.groupby(["seqname", f"peak_location_group{tolerance_type}"])["peak_location_actual"].transform('mean').astype(int)
        
        

      #   # set peak locations within LeftAlign_bp of eachother to eachother
      #   LeftAlign_bp = 80
      #   

      #   for i in range(len(peaks["file"].unique())):
            
      #       peaks["peak_location_ranges_diff_by_rep"] =peaks.groupby("seqname")["peak_location_shifted"].diff()
      #       peaks['peak_location_ranges_diff_by_rep'].fillna(0, inplace=True)
      #       peaks['peak_location_ranges_diff_by_rep'].replace(to_replace=0, method='ffill', inplace=True)

      #       peaks["peak_location_LeftAlign"] = peaks.apply(
      #          lambda row : row["peak_location_shifted"]-row["peak_location_ranges_diff_by_rep"] \
      #           if abs(row["peak_location_ranges_diff_by_rep"])<= self.peak_shift_value(row, self.peak_shift_value(row, LeftAlign_bp)) \
      #            else row["peak_location_shifted"], 
      #             axis =1
      #                      )
      #       peaks['peak_location_Centered'] = (
      #          peaks.groupby(["seqname", "peak_location_LeftAlign"])["peak_location_shifted"]
      #          .transform('mean')  # Calculate the mean
      #          .apply(lambda x: round(x / 10) * 10)  # Round to the nearest 10th
      #          .astype(int)  # Ensure the result is integer
      #       )
      #       peaks["peak_location_shifted"]=peaks["peak_location_Centered"] 
         
        
      #   peaks["peak_location"]=peaks.apply(lambda row:\
      #       row["peak_location_Centered"] \
      #          if abs(row["peak_location_actual"]- row["peak_location_Centered"]) < self.peak_shift_value(row, LeftAlign_bp)\
      #             else row["peak_location_actual"],
      #       axis=1)
        
      #   # one last diff to get the differences from the right
      #   peaks.sort_values(["seqname", "peak_location"], ascending=False, inplace=True)
      #   peaks["peak_location_ranges_diff_by_rep_right"] = abs(peaks.groupby("seqname")["peak_location_shifted"].diff())
      #   peaks['peak_location_ranges_diff_by_rep_right'].fillna(0, inplace=True)
      #   peaks['peak_location_ranges_diff_by_rep_right'].replace(to_replace=0, method='bfill', inplace=True)
      #   peaks.sort_values(["seqname", "peak_location"],ascending=True, inplace=True)
        
        
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

        peaks["peak_location_counts"] = peaks.groupby(["seqname",  "peak_location"])["peak_location"].transform('count')
        

        return peaks
    
     def finalize_call_peaks(self, peaks):
        """
        Applies functions to apply thresholds to data


        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files above a thresholds,
        
        """
        
        peaks = self.shift_peak_location(peaks)
        peaks = self.count_peak_occurances(peaks)

        return peaks

     def call_peaks(self):
        """
        Applies first derivative to bedgraph files to calculate the occurance of a peak based on passing through y=0
        by taking into account the location of neig peaks in the 1st derivative.
         determines order of max and min in derivative 
            peak is consecutive max and min
         |max|na| <-- derivative starts going up
         |na|min| <-- derivative ends going down
         # disregard noise peaks
         |na|min|   or   |max|na|
         |na|min|        |max|na|

        
        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files above a certain count threshold,
        
        """
        
        if os.path.isfile(self.WOnder_unfiltered_peaks):
            print("Reading in unfiltered peaks")
            peaks = pd.read_csv(self.WOnder_unfiltered_peaks)
            self.df_summary = pd.read_csv(os.path.join(self.WP_directory, "bedgraph_summary.csv"))
            return peaks

        print("Running WonderPeaks for ChIP...")
        pd.options.mode.chained_assignment = None
        
        # run funtions to call peaks
        bedfiles_dict = self.load_files()
        df = self.load_data(bedfiles_dict)
        df_peaks, df = self.first_derivative(df)
        

        df_peaks["next_min"] = df_peaks.groupby(["file","seqname"])["min"].shift(periods =-1)
        df_peaks["previous_max"] = df_peaks.groupby(["file","seqname"])["max"].shift(periods =1)
        df_peaks["previous_stop"] = df_peaks.groupby(["file","seqname"])["stop"].shift(periods =1)
        


        df_peaks.reset_index(inplace=True, drop =True)
        df_peaks_compared = df_peaks[
            ((~df_peaks["next_min"].isna()) & (~df_peaks["max"].isna()))
            |
            ((~df_peaks["previous_max"].isna()) & (~df_peaks["min"].isna()))
            |
            ((~df_peaks["previous_max"].isna()) & (~df_peaks["min"].isna()) & (~df_peaks["next_min"].isna()) & (~df_peaks["max"].isna()))
            ]



        df_peaks_compared.reset_index(inplace=True, drop =True)

        df_peaks_compared["peak_diff"] = df_peaks_compared.groupby(["file","seqname"])["start"].diff()
        df_peaks_compared.dropna(subset = ["min"], inplace=True)
        df_peaks_compared["peak_location"] = df_peaks_compared["stop"] - df_peaks_compared["peak_diff"] + 2*self.n
        df_peaks_compared["left_edge"] = df_peaks_compared["peak_location"] - df_peaks_compared["peak_diff"] - 2*self.n
        df_peaks_compared["right_edge"] = df_peaks_compared["peak_location"] + df_peaks_compared["peak_diff"] + 2*self.n
        
        
        peaks = pd.merge(df_peaks_compared, df, left_on =["seqname", "file", "peak_location"], 
                                                right_on = ["seqname","file","start"], how = "inner", suffixes  =("", "_peak"))
        peaks = pd.merge(peaks, df[["seqname","file","start", "score"]], left_on =["seqname", "file", "left_edge"], 
                                                right_on = ["seqname","file","start"], how = "inner", suffixes  =("", "_left_edge"))
        peaks = pd.merge(peaks, df[["seqname","file","start", "score"]], left_on =["seqname", "file", "right_edge"], 
                                                right_on = ["seqname","file","start"], how = "inner", suffixes  =("", "_right_edge"))
      #   peaks["score_left_edge"] = peaks["score_left_edge"].fillna(1)
      #   peaks["right_edge"] = peaks["right_edge"].fillna(1)
      #   peaks["Enrichment_left"] = peaks["score_peak"]/max(1, peaks["score_left_edge"])
      #   peaks["Enrichment_right"] = peaks["score_peak"]/max(1, peaks["score_right_edge"])
        
        
      #   peaks = pd.merge(self.df_summary, df_peaks_compared_merged, on = ["file", "seqname"])

        peaks = peaks[peaks["score"]>peaks["q75"]].copy()

                        
      #   peaks = self.finalize_call_peaks(peaks = peaks)
        peaks = self.shift_peak_location(peaks)
        peaks = self.count_peak_occurances(peaks)
        peaks.to_csv(self.WOnder_unfiltered_peaks, index = False)
        return peaks