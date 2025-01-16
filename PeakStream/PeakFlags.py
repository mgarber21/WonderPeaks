import os
import sys

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


# from Helper_Functions import *
from PeakStream.AT_stretches import *
from PeakStream.WonderPeaks4UTRs import *
from PeakStream.PeakLinks import *

# Parent_directory = os.path.dirname(os.getcwd())  # Parent file's directory
# sys.path.insert(0, Parent_directory)
from .Helper_Functions import *


class PeakStreamFlag(PeakStreamLinks):
    
    def __init__(self,  directory, coordinate_file, 
                 genome_fasta_dir,  n, stretch_length):
        
        ATstretches.__init__(self,directory, coordinate_file, genome_fasta_dir, stretch_length)
        PeakStreamLinks.__init__(self, directory, coordinate_file, 
                                 genome_fasta_dir,  n, stretch_length)

        """
        Initialize cleaning peaks links.

        Parameters:
        - directory (str): path to directory with data
        - coordinate_file (str): path to genome coordinate file
        - n (int): steps used to calculate peaks
        - genome_fasta_dir (str): path to genome fasta file
        - stretch_length (int): cutoff length of A/T stretches

        """

        self.AT_coordinates_inGene = self.ApplyATStretches()
        self.peak_links = self.ApplyPeakLinks()
        self.df_summary = pd.read_csv(os.path.join(self.PS_directory, "bedgraph_summary.csv"))
        self.AT_Flags_output = os.path.join(self.PS_directory, "Peaks2ATstretches.csv")
        self.AT2peaks_output = os.path.join(self.PS_directory, "ATstretches2Peaks.csv")
        


    def Link_w_AT(self):
        """
        link peaks to nearest AT-stretches within genes
        
        return:
        -pd.DataFrame: DataFrame of coordinates of A/T stretches matched to nearest peak

        """

        #use output from this to indicate that peaks are not reliable poly-A-primed events
        self.peak_links['peak_location'] = self.peak_links['peak_location'].astype("int")
        self.peak_links["nucleotide"] = self.peak_links.apply(lambda row: "A" if row["strand"] == "+" else "T", axis=1)
        self.AT_coordinates_inGene['first'] = self.AT_coordinates_inGene['first'].astype("int")
        self.peak_links.sort_values(by ="peak_location", inplace=True)
        self.AT_coordinates_inGene.sort_values(by ="first", inplace=True)

        tolerance = 200
        self.peak_links["index"] = self.peak_links.index

        peak_links_AT = pd.DataFrame()
        for method, strand in zip(["forward", "backward"], ["+", "-"]):
            df_linked_peaks_nucleotide = self.peak_links[self.peak_links["strand"] == strand]

            # here method direction is linked to nucleotide and strand
            peak_links_AT_temp = pd.merge_asof( df_linked_peaks_nucleotide, self.AT_coordinates_inGene, 

                        left_index=False, right_index=False, 
                        left_on="peak_location",  right_on="first", 
                        by=["seqname", "strand" ],
                        suffixes=('','_AT'), tolerance=tolerance, 
                        allow_exact_matches=True, direction=method)
            peak_links_AT = pd.concat([peak_links_AT, peak_links_AT_temp])

        peak_links_AT.drop_duplicates(subset =[ "seqname", "peak_location", "first", "file"], inplace=True)
        peak_links_AT.dropna(how ="any", subset =["seqname","peak_location", "first"], inplace=True)
        
        peak_links_AT.rename(
            columns = dict(zip([col for col in peak_links_AT.columns\
                                if col not in self.peak_links.columns and "AT" not in col], 
                                [col+"_AT" for col in peak_links_AT.columns\
                                if col not in self.peak_links.columns and "AT" not in col])),
            inplace=True)

        return peak_links_AT

    def AT_Stretch(self, peak_links_AT):
        """
        link peaks linked to genes to nearest AT-stretches within genes
        
        return:
        -pd.DataFrame: peak_links DataFrame with bool column indicating association with A/T stretch
        
        """
        print("Assigning peaks to A/T stretches within gene boundaries...")
        peak_links = self.peak_links



        # whether a peak is truly associated with an A/T stretch, probably depends on the 
        # gene size and the location of the peak relative to the end of the gene
        # - gene_size: length of the gene from start to end
        # - gene_factor: the smaller of 500 or half the gene size
        # - AT_strech: whether the peaks meets the follwing conditions:
        #   1. the start of the gene linked to the A/T stretch matches that of the peak
        #   2. If the strand is positive or negative, the peak lies between the start of the gene and near the end of the gene
        #   2a. how close the peak can be to the end of the gene depends on the gene_factor
        peak_links_AT["gene_size"] = abs(peak_links_AT["start"]-peak_links_AT["end"]).astype(int)
        peak_links_AT["gene_factor"] = np.minimum(500, peak_links_AT["gene_size"]*0.5).astype(int)


        peak_links_AT["AT_stretch"] = (
                        (peak_links_AT["start_AT"] == peak_links_AT["start"])  
                        &
                        (
                            ((peak_links_AT["strand"] == "+")
                            &
                            (peak_links_AT["peak_location"].between(
                                peak_links_AT["start_AT"] - 50,peak_links_AT["end_AT"] - peak_links_AT["gene_factor"]
                                )
                                ))
                                |
                                ((peak_links_AT["strand"] == "-")
                                &
                                (peak_links_AT["peak_location"].between(
                                    peak_links_AT["start_AT"] + peak_links_AT["gene_factor"] ,peak_links_AT["end_AT"] + 50
                                    )
                                    ))))
        # save At-->peaks file for analysis
        peak_links_AT.to_csv(self.AT2peaks_output, index=False)

        # which peaks should not be considered because of their A/T stretch association
        unlinked_peaks = peak_links_AT[peak_links_AT["AT_stretch"]]["index"]


        # which peaks are identified as having AT-stretch association, and append info the the peaks dataframe
        AT_stretch = (peak_links["index"].isin(list(unlinked_peaks)))
        self.peak_links["AT_Stretch"] = AT_stretch
        return self.peak_links


        

    def blacklist(self, peak_links_AT):
        """
        blacklist genes that have highscoring peaks in range of AT_rich region in a tandem gene
        these genes should be excluded from determining peak boundaries
        
        return:
        -pd.DataFrame: peak_links DataFrame with bool column indicating association with blacklist
        
        """
        
        print("Assigning blacklisted (bl) peaks...")
        peak_links = self.peak_links
        peak_links_AT = pd.merge(self.df_summary, peak_links_AT, on = "file")

        # blacklist genes that have highscoring peaks in range of AT_rich region
        # these genes should be excluded from determining peak boundaries
        
        # identify peaks downstream of AT-stretches that fall within the next tandem gene
        black_list = ((peak_links_AT["score_peak"] > peak_links_AT["q90"]) 
                    & (((peak_links_AT["start"] == peak_links_AT["start_preceding_gene_AT"]) & (peak_links_AT["strand"] == "+") )
                    | ((peak_links_AT["start"] == peak_links_AT["start_next_gene_AT"]) & (peak_links_AT["strand"] == "-") ))
                    
                    & (peak_links_AT.duplicated(subset = ["file", "seqname", "strand", "start"], keep =False)))
        

        # create a dataframe with linking blacklisted peaks to AT-stretch dataframe 
        peak_links_AT_blacklist = pd.DataFrame()
        for strand, next_gene_start in zip(["+", "-"], ["start_preceding_gene_AT", "start_next_gene_AT"]):
            blacklist_temp = peak_links_AT[(black_list) & (peak_links_AT["strand"] == strand)][["file", "seqname", "strand",
                                                                                                                            "start_AT", "end_AT",
                                                                                                                            "start_preceding_gene_AT", "end_preceding_gene_AT", 
                                                                                                                            "start_next_gene_AT", "end_next_gene_AT",
                                                                                                                            "peak_location"]]
            blacklist_temp["start"] = blacklist_temp[next_gene_start]
            blacklist_temp["blacklist"] = True
            peak_links_AT_blacklist = pd.concat([blacklist_temp, peak_links_AT_blacklist])


        # merge the blacklist-AT dataframe with the comeplete peaks dataframe
        peak_links["strand_sign"] = peak_links.apply(lambda row: 1 if row["strand"]== "+" else -1, axis =1)
        peak_links_blacklist = pd.merge(peak_links, peak_links_AT_blacklist,  on = ["file", "seqname", "strand", "start"], suffixes=("", "_blacklist"), how = "outer")
        peak_links_blacklist.fillna({"blacklist":False, "peak_location_blacklist":peak_links_blacklist["peak_location"] - peak_links_blacklist["strand_sign"]*100 }, inplace=True)

        # identify blacklisted peaks that fall downstream of AT-stretch associated peak
        peak_links_blacklist_pos = peak_links_blacklist[
            ((peak_links_blacklist["blacklist"]) & ((peak_links_blacklist["strand"] == "+") & (peak_links_blacklist["peak_location"] >=  peak_links_blacklist["peak_location_blacklist"])))
            ]
        peak_links_blacklist_neg = peak_links_blacklist[
            ((peak_links_blacklist["blacklist"]) & ((peak_links_blacklist["strand"] == "-") & (peak_links_blacklist["peak_location"] <= peak_links_blacklist["peak_location_blacklist"])))
            ]
        peak_links_blacklist = pd.concat([peak_links_blacklist_pos, peak_links_blacklist_neg])

        # merge the blacklist dataframe with the comeplete peaks dataframe
        self.peak_links = pd.merge(peak_links, peak_links_blacklist[["file", "seqname", "strand", "start","blacklist",  "peak_location", "peak_location_blacklist"]],
                on = ["seqname", "strand", "start", "peak_location"], suffixes=("", "_bl"), how = "left").drop_duplicates(["file", "seqname", "strand", "start","blacklist",  "peak_location"], keep = "first")
        self.peak_links.fillna({"blacklist":False, "peak_location_blacklist":peak_links["peak_location"] - peak_links["strand_sign"]*100 }, inplace=True)

        
    
    def ApplyFlags(self):

        """
        run codes within PeakStreamFlag

        return"
        -pd.DataFrame: peak_links datafram with A/T strech and blacklist flags
        
        """
        

        peak_links_AT = self.Link_w_AT()
        self.AT_Stretch(peak_links_AT)
        self.blacklist(peak_links_AT)

        
        self.peak_links.to_csv(self.AT_Flags_output, index=False)


        return self.peak_links