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

# Parent_directory = os.path.dirname(os.getcwd())  # Parent file's directory
# sys.path.insert(0, Parent_directory)
from .Helper_Functions import *


class PeakStreamLinks(ATstretches, WonderPeaks_4UTRs):

    """
    Link peaks to gene coordinates    
    
    """
    
    def __init__(self, directory, coordinate_file, 
                  genome_fasta_dir,  n, stretch_length):
        ATstretches.__init__(self,directory, coordinate_file, genome_fasta_dir, stretch_length)
        WonderPeaks_4UTRs.__init__(self, directory, coordinate_file)

        """
        Initialize linking peaks to genes.

        Parameters:
        - directory (str): path to directory with data
        - coordinate_file (str): path to genome coordinate file
        - n (int): steps used to calculate peaks
        - genome_fasta_dir (str): path to genome fasta file
        - stretch_length (int): cutoff length of A/T stretches
        """

        self.peaks = self.call_peaks()
        self.AT_coordinates_inGene = self.ApplyATStretches()
        
    def link_peaks_to_GTF(self):
        print("Linking peaks to gene coordinates...")
        
        self.peaks ['peak_location'] = self.peaks ['peak_location'].astype("int")
        self.peaks.sort_values(by ="peak_location", inplace=True)

        self.gtf_data['start'] = self.gtf_data['start'].astype("int")
        self.gtf_data.sort_values(by ="start", inplace=True)

        cols_peaks = ['seqname', 'start', 'stop', 'score', 'file',
                    'strand', 'nucleotide', 'peak_location', 'start_peak', 'stop_peak', 'score_peak',
                    'strand_direction', "peak_diff"]


        peaks_GTF = pd.DataFrame()
        
        for gene_position in ["start", "end"]:
            for method in ["forward", "backward"]:

                self.gtf_data.sort_values(by =gene_position, inplace=True)
                peaks_GTF_temp = pd.merge_asof(self.peaks[cols_peaks], self.gtf_data[self.cols],
                            # on="seqname", 

                            left_index=False, right_index=False, 
                            left_on="peak_location", right_on=gene_position, 
                            by=["seqname", "strand" ],
                            suffixes=('_data', ''), tolerance=20000, 
                            allow_exact_matches=True, direction=method)
                peaks_GTF_temp["method"] = method
                peaks_GTF_temp["on"] = gene_position
                peaks_GTF = pd.concat([peaks_GTF, peaks_GTF_temp])

        # drop anything that is unlinked
        peaks_GTF.dropna(how ="any", subset =["seqname","peak_location", "start", "end"], inplace=True)

        # negative must be linked to start of genes and positve must be linked to end of genes
        neg_method = (peaks_GTF["on"] == "start") & (peaks_GTF["strand"] == "-")
        pos_method = (peaks_GTF["on"] == "end") & (peaks_GTF["strand"] == "+") 
        peaks_GTF = peaks_GTF[neg_method | pos_method]


        ## some additional data cleanup for PeakStream
        for gene_position in ["start", "end"]:
            peaks_GTF[f"peak_to_{gene_position}"] = peaks_GTF["peak_location"] - peaks_GTF[gene_position]


        peaks_GTF.sort_values(["file","seqname", "peak_location","start", "end", ], inplace=True)
        peaks_GTF.fillna(dict(zip(['end_next_gene', 'end_preceding_gene', 'start_next_gene', 'start_preceding_gene'], [1]*4)), inplace=True)
        peaks_GTF["gene_length"] = abs(peaks_GTF["start"] - peaks_GTF["end"])
        peaks_GTF["next_gene_length"] = abs(peaks_GTF["start_next_gene"] - peaks_GTF["end_next_gene"])
        peaks_GTF["preceding_gene_length"] = abs(peaks_GTF["start_preceding_gene"] - peaks_GTF["end_preceding_gene"])

        # linked peaks are not allowed to be downstream of next gene
        peaks_GTF =  peaks_GTF[(
                                (peaks_GTF["peak_location"].between(peaks_GTF['start'], peaks_GTF['start_next_gene']))
                                &
                                (peaks_GTF["strand"] == "+")
                                )
                                |
                                (
                                    (peaks_GTF["peak_location"].between(peaks_GTF['end_preceding_gene'],peaks_GTF['end']))
                                    &
                                    (peaks_GTF["strand"] == "-")
                                    )
                                    ]
                        
        
        
        return peaks_GTF

    def file_stream(self):
        """
        Apply stream function to link peaks to genes for each file

        returns: 
        - dict: dictionary of bedgraph file to dictionary of peaks linked to genes
        """
        peaks_GTF = self.link_peaks_to_GTF()
        file_peak_links_dict = dict()
        for file in peaks_GTF["file"].unique():
            peaks_GTF_file = peaks_GTF[peaks_GTF["file"] == file]

            peak_links_dict = dict()
            

            for i, row in peaks_GTF_file.iterrows():

                dict_peak_link = stream(row)
                peak_links_dict.update(dict_peak_link)
        
        
            file_peak_links_dict[file] = pd.DataFrame.from_dict(peak_links_dict, orient='index', columns= ["seqname", "strand", "start", "end", "score_peak"]).reset_index().rename(columns={"index":"peak_location"})

        return file_peak_links_dict

    def ApplyPeakLinks(self):
        """
        Apply file_stream function and create a datafram from the file-->data dictionary

        returns: 
        - pd.DataFrame: DataFrame of peaks linked to genes for each bedgraph file
        """
        
        file_peak_links_dict = self.file_stream()
        peak_links = pd.concat(file_peak_links_dict).droplevel(1).reset_index().rename(columns={"index":"file"})

        return peak_links
