"""
PeakStream is a Python-based approach to defining transcript boundaries
in poly-A-priming library RNA-seq data.

The approach uses WonderPeaks to call peaks based on the first derivative
and links them to known transcripts or identifies hypothetical transcripts.
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


# 25th Percentile
def q25(x):
    return x.quantile(0.25)
# 50th Percentile
def q50(x):
    return x.quantile(0.5)
# 75th Percentile
def q75(x):
    return x.quantile(0.75)
def q90(x):
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

    dfGTF_GeneRef = pd.read_csv(coordinate_file, index_col=False)
    dfGTF_GeneRef.sort_values(["seqname", "start","end"], inplace = True)

    dfGTF_GeneRef["start_next_gene"] = dfGTF_GeneRef.groupby(["seqname", "strand"])["start"].shift(periods =-1)
    dfGTF_GeneRef["strand_next_gene"] = dfGTF_GeneRef.groupby(["seqname", "strand"])["strand"].shift(periods =-1)
    dfGTF_GeneRef["end_next_gene"] = dfGTF_GeneRef.groupby(["seqname", "strand"])["end"].shift(periods =-1)


    dfGTF_GeneRef["end_preceding_gene"] = dfGTF_GeneRef.groupby(["seqname", "strand"])["end"].shift(periods =1)
    dfGTF_GeneRef["strand_preceding_gene"] = dfGTF_GeneRef.groupby(["seqname", "strand"])["strand"].shift(periods =1)
    dfGTF_GeneRef["start_preceding_gene"] = dfGTF_GeneRef.groupby(["seqname", "strand"])["start"].shift(periods =1)

    return dfGTF_GeneRef

class AT_stretches():
    def __init__(self, coordinate_file, genome_fasta_dir, stretch_length):
        """
        Initialize and load genome sequences and annotations for A/T stretches.

        Parameters:
        - coordinate_file (str): path to genome coordinate file
        - genome_fasta_dir (str): path to genome fasta file
        - stretch_length (int): cutoff length of A/T stretches
        """

        self.coordinate_file  = coordinate_file
        self.genome_fasta_dir = genome_fasta_dir 
        self.stretch_length = stretch_length
        self.next_gene_cols = ["start_next_gene", "strand_next_gene","end_next_gene",
                            "end_preceding_gene","strand_preceding_gene","start_preceding_gene"]
        self.cols = ['seqname', 'start', 'end',  'strand']+self.next_gene_cols
        self.dfGTF_GeneRef  = load_coordinates(coordinate_file = coordinate_file)


    def find_AT_stretches(self):
         """
        Identify stretches of A/T nucleotides exceeding a defined length.
        """

        with open(self.genome_fasta_dir) as fasta_file:  # Will close handle cleanly
            identifiers = []
            seqs = []
            iter = []
            for title, sequence in SimpleFastaParser(fasta_file):
                identifiers.append(title.split(None, 1)[0])  # First word is ID
                seqs.append(list(sequence))
                iter.append(list(range(len(sequence))))

        df_genome_seqs = pd.DataFrame([identifiers, seqs, iter]).T.explode([1, 2]).reset_index()
        df_genome_seqs.rename(columns={0:"seqname", 1:"nucleotide", 2:"iter"}, inplace=True)
        df_genome_seqs['subgroup'] = (df_genome_seqs['nucleotide'] != df_genome_seqs['nucleotide'].shift(1)).cumsum()
        df_genome_seqs_duplicated = df_genome_seqs[df_genome_seqs["subgroup"].duplicated(keep = False)]
        subgroup_counts = df_genome_seqs_duplicated['subgroup'].value_counts()
        df_genome_seqs_duplicated_5X = df_genome_seqs_duplicated[df_genome_seqs_duplicated['subgroup'].isin(subgroup_counts[subgroup_counts>=self.stretch_length].index)]
        df_genome_seqs_duplicated_5X_grouped = df_genome_seqs_duplicated_5X.groupby(["seqname", "nucleotide", "subgroup"],as_index=False).agg({"iter":["first", "last","count"]}).droplevel(axis=1, level = 0)
        df_genome_seqs_duplicated_5X_grouped.columns = ["seqname", "nucleotide", "subgroup"] +['first', 'last', 'count']
        df_genome_seqs_duplicated_5X_grouped_AT= df_genome_seqs_duplicated_5X_grouped[df_genome_seqs_duplicated_5X_grouped["nucleotide"].isin(["A", "T"])]
    
        return df_genome_seqs_duplicated_5X_grouped_AT
    
    def plot_AT_dist(self):
        df_genome_seqs_duplicated_5X_grouped_AT = self.find_AT_stretches()
        fig, ax = plt.subplots(figsize = (15,5))
        sns.barplot(data = pd.DataFrame(df_genome_seqs_duplicated_5X_grouped_AT.value_counts(["count", "nucleotide"])).reset_index(),
                    x="count",y=0,
                    hue = "nucleotide", ax = ax, palette= ["#000000", "#000587"],
                    linewidth=2.5, edgecolor="#ffffff"
                    )
        ax.set_yscale("log")
        ax.set_ylabel("# of A ot T stretches of {x} length")
        ax.set_xlabel("Length of A or T stretches")

        ax.bar_label(ax.containers[0], fontsize=10, padding = 1, color ="#000000" )
        ax.bar_label(ax.containers[1], fontsize=10, padding = 10, color ="#000587")
        ax.set_xlim(-1, 27)
        ax.set_ylim(ax.get_ylim()[0], 10**5)
        ax.set_title("Distribution of A or T stretches across the genome")
        ax.legend(loc = "upper right", ncol =2, frameon =False)
        sns.despine(top =True, right = True)
        plt.show()
    
    def merge_AT_w_GTF(self):

        df_genome_seqs_duplicated_5X_grouped_AT = self.find_AT_stretches()


        df_genome_seqs_duplicated_5X_grouped_AT['first'] = df_genome_seqs_duplicated_5X_grouped_AT['first'].astype("int")
        self.dfGTF_GeneRef['start'] = self.dfGTF_GeneRef['start'].astype("int")
        df_genome_seqs_duplicated_5X_grouped_AT.sort_values(by ="first", inplace=True)
        self.dfGTF_GeneRef.sort_values(by ="start", inplace=True)
        df_genome_seqs_duplicated_5X_grouped_AT["strand"] = df_genome_seqs_duplicated_5X_grouped_AT.apply(lambda row: "+" if row["nucleotide"] == "A" else "-", axis= 1)


        df_genome_seqs_duplicated_5X_grouped_AT_GTF = pd.DataFrame()


        for method in ["backward", "forward"]:
            self.dfGTF_GeneRef["start"] = self.dfGTF_GeneRef["start"].astype("int")
            self.dfGTF_GeneRef.sort_values(by ="start", inplace=True)
            df_genome_seqs_duplicated_5X_grouped_AT_GTF_temp = pd.merge_asof(df_genome_seqs_duplicated_5X_grouped_AT, self.dfGTF_GeneRef[self.cols], 
                        left_index=False, right_index=False, 
                        left_on=["first"], right_on="start", 
                        by=["seqname", "strand"],
                        suffixes=('_x', '_y'), tolerance=20000, 
                        allow_exact_matches=True, direction=method)
            df_genome_seqs_duplicated_5X_grouped_AT_GTF = pd.concat([df_genome_seqs_duplicated_5X_grouped_AT_GTF, df_genome_seqs_duplicated_5X_grouped_AT_GTF_temp])

        df_genome_seqs_duplicated_5X_grouped_AT_GTF.dropna(how ="any", inplace=True)
        df_genome_seqs_duplicated_5X_grouped_AT_GTF.drop_duplicates(subset =["seqname", "first", "last", "start", "end"], inplace=True)

        return df_genome_seqs_duplicated_5X_grouped_AT_GTF
    

    def reduce_AT_GTF(self):

        df_genome_seqs_duplicated_5X_grouped_AT_GTF = self.merge_AT_w_GTF()
        # AT regions witin the first 85% of a gene if the gene >

        #gene:             ------>
        #ATcoverage:       -----
        pos_AT = ((df_genome_seqs_duplicated_5X_grouped_AT_GTF["strand"] == "+"))
        neg_AT = ((df_genome_seqs_duplicated_5X_grouped_AT_GTF["strand"] == "-"))

        df_genome_seqs_duplicated_5X_grouped_AT_GTF_within =df_genome_seqs_duplicated_5X_grouped_AT_GTF[pos_AT|neg_AT]
        within_gene = (
            (df_genome_seqs_duplicated_5X_grouped_AT_GTF_within["first"] > df_genome_seqs_duplicated_5X_grouped_AT_GTF_within["start"])
                & (df_genome_seqs_duplicated_5X_grouped_AT_GTF_within["first"] < df_genome_seqs_duplicated_5X_grouped_AT_GTF_within["end"])
                )
        df_genome_seqs_duplicated_5X_grouped_AT_GTF_within = df_genome_seqs_duplicated_5X_grouped_AT_GTF_within[within_gene]

        return df_genome_seqs_duplicated_5X_grouped_AT_GTF_within
    

class WonderPeaks_UTRs():
    """
        Process bedgraph files and detect transcript peaks using the first derivative.
        """

    def __init__(self, coordinate_file, bed_directory, n):

        """
        Initialize and load genome sequences and bedgraph files for peak calling.

        Parameters:
        - coordinate_file (str): path to genome coordinate file
        - bed_directory (str): path to directory with bedgraph files
        - n (int): steps used to calculate peaks
        """


        self.dfGTF_GeneRef  = load_coordinates(coordinate_file)
        self.n = n
        self.bed_directory = bed_directory

        

    def load_files(self):
        
        # add all bedgraph files relevant to the study to a single folder

        strands = ["fwd", "rev"]
        files = []
        dict_files  = {}
        for strand in strands:
            
            files = []
            for file in [file for file in os.listdir(self.bed_directory) if strand in file]:
                files.append( os.path.join(self.bed_directory,file))
                files = sorted(files)
            dict_files[strand] = files

        return dict_files
    
    def load_data(self):
        dict_files = self.load_files()

        global df
        df = pd.DataFrame()
        for strand, files in dict_files.items():
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

        global df_summary
        
        df_summary = df.groupby("file")["score"].agg(["median","mean", q25, q50, q75, q90])


        return df
    
    def first_der(self):

        # n =10  # number of points to be checked before and after the peak. Using 10 bc its half the binsize (20).    
        
        df = self.load_data()

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
    

    def median_score_cutoff(self, df_peaks_all_comparisons_merge):

        df_peaks_all_comparisons_merge_score_cutoff = pd.merge(df_peaks_all_comparisons_merge, df_summary, left_on = "file", right_index=True)

        df_peaks_all_comparisons_merge_score_cutoff= df_peaks_all_comparisons_merge_score_cutoff[
            df_peaks_all_comparisons_merge_score_cutoff["score_peak"] > df_peaks_all_comparisons_merge_score_cutoff["median"]
            ]

        return  df_peaks_all_comparisons_merge_score_cutoff

    def clean_call_peaks(self, df_all_peaks):

        df_all_peaks['nucleotide'] = df_all_peaks.apply(lambda row: "A" if row['strand'] == "fwd" else "T", axis= 1)
        df_all_peaks['dir'] = df_all_peaks['strand']
        df_all_peaks['strand'] = df_all_peaks.apply(lambda row: "+" if row['strand'] == "fwd" else "-", axis= 1)
        df_all_peaks.rename(columns={"chr":"seqname"}, inplace =True)
        df_all_peaks[df_all_peaks["peak_diff"] < 500]

        return df_all_peaks
    


    def normalize_peak_location(self, df_all_peaks):

        # set peak locations within 30bp of eachother to eachother

        df_all_peaks.sort_values(["seqname", "peak_location"], inplace=True)
        df_all_peaks["peak_location_ranges_diff_by_rep"] =df_all_peaks.groupby("seqname")["peak_location"].diff()
        df_all_peaks['peak_location_ranges_diff_by_rep'].fillna(0, inplace=True)

        df_all_peaks['peak_location_ranges_diff_by_rep'].replace(to_replace=0, method='ffill', inplace=True)

        df_all_peaks["peak_location_LeftAlign"] = df_all_peaks.apply(lambda row : row["peak_location"]-row["peak_location_ranges_diff_by_rep"] \
                                                                     if abs(row["peak_location_ranges_diff_by_rep"])<= 100 \
                                                                        else row["peak_location"], 
                                                                        axis =1)
        df_all_peaks["peak_location_Centered"] = df_all_peaks.groupby(["seqname","strand",  "peak_location_LeftAlign"])["peak_location"].transform('mean')

        df_all_peaks["peak_location"]=df_all_peaks["peak_location_Centered"] 

        return df_all_peaks

    def count_peak_occurances(self, df_all_peaks):

        # count the number of occurances of a peak_location
        # peak must occur more than twice to be considered for linking

        # peak_location_counts = pd.DataFrame(df_all_peaks[["seqname","strand", "peak_location"]].value_counts()).rename(columns = {0:"peak_location_counts"}).reset_index()
        df_all_peaks["peak_location_counts"] = df_all_peaks.groupby(["seqname","strand",  "peak_location"])["peak_location"].transform('count')
        # df_all_peaks = pd.merge(df_all_peaks, peak_location_counts, on =["seqname", "strand", "peak_location"])
        count_centeredpeaks = ((df_all_peaks["peak_location_counts"] > 2))
        df_all_peaks = df_all_peaks[count_centeredpeaks]

        return df_all_peaks
    
    def finalize_call_peaks(self, df_all_peaks):
        df_all_peaks = self.clean_call_peaks(df_all_peaks)
        df_all_peaks = self.normalize_peak_location(df_all_peaks)
        df_all_peaks = self.count_peak_occurances(df_all_peaks)

        return df_all_peaks

    def call_peaks(self, median_score_cutoff_ = True):
        pd.options.mode.chained_assignment = None

        df_peaks, df = self.first_der()
        

        df_peaks["next_min"] = df_peaks.groupby(["file","chr"])["min"].shift(periods =-1)
        df_peaks["previous_max"] = df_peaks.groupby(["file","chr"])["max"].shift(periods =1)
        df_peaks["previous_stop"] = df_peaks.groupby(["file","chr"])["stop"].shift(periods =1)
        df_peaks.reset_index(inplace=True, drop =True)
        df_peaks_all_comparisons = df_peaks[
            ((~df_peaks["next_min"].isna()) & (~df_peaks["max"].isna()))
            |
            ((~df_peaks["previous_max"].isna()) & (~df_peaks["min"].isna()))
            |
            ((~df_peaks["previous_max"].isna()) & (~df_peaks["min"].isna()) & (~df_peaks["next_min"].isna()) & (~df_peaks["max"].isna()))
            ]



        df_peaks_all_comparisons.reset_index(inplace=True, drop =True)

        df_peaks_all_comparisons["peak_diff"] = df_peaks_all_comparisons.groupby(["file","chr"])["start"].diff()
        df_peaks_all_comparisons.dropna(subset = ["min"], inplace=True)
        df_peaks_all_comparisons["peak_location"] = df_peaks_all_comparisons["stop"] - df_peaks_all_comparisons["peak_diff"] + 2*self.n
        df_peaks_all_comparisons_merge = pd.merge(df_peaks_all_comparisons, df, left_on =["chr", "file", "peak_location"], 
                                                right_on = ["chr","file","start"], how = "inner", suffixes  =("", "_peak"))

        if median_score_cutoff_:
            df_peaks_all_comparisons_merge_score_cutoff = self.median_score_cutoff(df_peaks_all_comparisons_merge)
            return self.finalize_call_peaks(df_all_peaks = df_peaks_all_comparisons_merge_score_cutoff)
        

        return self.finalize_call_peaks(df_all_peaks = df_peaks_all_comparisons_merge)
    

class PeakStreamLinks(AT_stretches, WonderPeaks_UTRs):
    
    def __init__(self, bed_directory, coordinate_file, 
                  genome_fasta_dir,  n, stretch_length):
        AT_stretches.__init__(self, coordinate_file, genome_fasta_dir, stretch_length)
        WonderPeaks_UTRs.__init__(self, coordinate_file, bed_directory, n)

                """
        Initialize linking peaks to genes.

        Parameters:
        - bed_directory (str): path to directory with bedgraph files
        - coordinate_file (str): path to genome coordinate file
        - n (int): steps used to calculate peaks
        - genome_fasta_dir (str): path to genome fasta file
        - stretch_length (int): cutoff length of A/T stretches
        """

        self.df_all_peaks = self.call_peaks()
        self.df_genome_seqs_duplicated_5X_grouped_AT_GTF_within = self.reduce_AT_GTF()
        
    def link_peaks_to_GTF(self):
        
        self.df_all_peaks ['peak_location'] = self.df_all_peaks ['peak_location'].astype("int")
        self.dfGTF_GeneRef['start'] = self.dfGTF_GeneRef['start'].astype("int")
        self.df_all_peaks.sort_values(by ="peak_location", inplace=True)
        self.dfGTF_GeneRef.sort_values(by ="start", inplace=True)

        cols_peaks = ['seqname', 'start', 'stop', 'score', 'file',
                    'strand', 'nucleotide', 'peak_location', 'start_peak', 'stop_peak', 'score_peak',
                    'dir', "peak_diff"]


        df_all_peaks_GTF = pd.DataFrame()
        for gene_position in ["start", "end"]:
            for method in ["forward", "backward"]:

                self.dfGTF_GeneRef.sort_values(by =gene_position, inplace=True)
                df_all_peaks_GTF_temp = pd.merge_asof(self.df_all_peaks[cols_peaks], self.dfGTF_GeneRef[self.cols],
                            # on="seqname", 

                            left_index=False, right_index=False, 
                            left_on="peak_location", right_on=gene_position, 
                            by=["seqname", "strand" ],
                            suffixes=('_data', ''), tolerance=20000, 
                            allow_exact_matches=True, direction=method)
                df_all_peaks_GTF_temp["method"] = method
                df_all_peaks_GTF_temp["on"] = gene_position
                df_all_peaks_GTF = pd.concat([df_all_peaks_GTF, df_all_peaks_GTF_temp])

        # drop anything that is unlinked
        df_all_peaks_GTF.dropna(how ="any", subset =["seqname","peak_location", "start", "end"], inplace=True)

        # negative must be linked to start of genes and positve must be linked to end of genes
        neg_method = (df_all_peaks_GTF["on"] == "start") & (df_all_peaks_GTF["strand"] == "-")
        pos_method = (df_all_peaks_GTF["on"] == "end") & (df_all_peaks_GTF["strand"] == "+") 
        df_all_peaks_GTF = df_all_peaks_GTF[neg_method | pos_method]


        ## some additional data cleanup for PeakStream
        for gene_position in ["start", "end"]:
            df_all_peaks_GTF[f"peak_to_{gene_position}"] = df_all_peaks_GTF["peak_location"] - df_all_peaks_GTF[gene_position]


        df_all_peaks_GTF.sort_values(["file","seqname", "peak_location","start", "end", ], inplace=True)
        df_all_peaks_GTF.fillna(dict(zip(['end_next_gene', 'end_preceding_gene', 'start_next_gene', 'start_preceding_gene'], [1]*4)), inplace=True)
        df_all_peaks_GTF["gene_length"] = abs(df_all_peaks_GTF["start"] - df_all_peaks_GTF["end"])
        df_all_peaks_GTF["next_gene_length"] = abs(df_all_peaks_GTF["start_next_gene"] - df_all_peaks_GTF["end_next_gene"])
        df_all_peaks_GTF["preceding_gene_length"] = abs(df_all_peaks_GTF["start_preceding_gene"] - df_all_peaks_GTF["end_preceding_gene"])

        return df_all_peaks_GTF

    def data_stream(self):
        df_all_peaks_GTF = self.link_peaks_to_GTF()
        file_peak_links_dict = dict()
        for file in df_all_peaks_GTF["file"].unique():
            df_all_peaks_GTF_file = df_all_peaks_GTF[df_all_peaks_GTF["file"] == file]

            peak_links_dict = dict()
            

            for i, row in df_all_peaks_GTF_file.iterrows():

                dict_peak_link = stream(row)
                peak_links_dict.update(dict_peak_link)
        
        
            file_peak_links_dict[file] = pd.DataFrame.from_dict(peak_links_dict, orient='index', columns= ["seqname", "strand", "start", "end", "score_peak"]).reset_index().rename(columns={"index":"peak_location"})

        return file_peak_links_dict

    def concat_data(self):
        file_peak_links_dict = self.data_stream()
        df_peak_links = pd.concat(file_peak_links_dict).droplevel(1).reset_index().rename(columns={"index":"file"})

        return df_peak_links
    
class PeakStreamClean(PeakStreamLinks):
    
    def __init__(self,  bed_directory, coordinate_file, 
                 genome_fasta_dir,  n, stretch_length):
        
        PeakStreamLinks.__init__(self, bed_directory, coordinate_file, 
                                 genome_fasta_dir,  n, stretch_length)

                """
        Initialize cleaning peaks links.

        Parameters:
        - bed_directory (str): path to directory with bedgraph files
        - coordinate_file (str): path to genome coordinate file
        - n (int): steps used to calculate peaks
        - genome_fasta_dir (str): path to genome fasta file
        - stretch_length (int): cutoff length of A/T stretches
        """

        self.df_genome_seqs_duplicated_5X_grouped_AT_GTF_within = self.reduce_AT_GTF()
        self.df_peak_links = self.concat_data()


    def Link_w_AT(self):
   
        # link all peaks to nearest AT-stretches within genes

        #use output from this to indicate that peaks are not reliable poly-A-primed events
        self.df_peak_links['peak_location'] = self.df_peak_links['peak_location'].astype("int")
        self.df_peak_links["nucleotide"] = self.df_peak_links.apply(lambda row: "A" if row["strand"] == "+" else "T", axis=1)
        self.df_genome_seqs_duplicated_5X_grouped_AT_GTF_within['first'] = self.df_genome_seqs_duplicated_5X_grouped_AT_GTF_within['first'].astype("int")
        self.df_peak_links.sort_values(by ="peak_location", inplace=True)
        self.df_genome_seqs_duplicated_5X_grouped_AT_GTF_within.sort_values(by ="first", inplace=True)

        tolerance = 200
        self.df_peak_links["index"] = self.df_peak_links.index

        df_linked_peaks_AT = pd.DataFrame()
        for method, strand in zip(["forward", "backward"], ["+", "-"]):
            df_linked_peaks_nucleotide = self.df_peak_links[self.df_peak_links["strand"] == strand]

            # here method direction is linked to nucleotide and strand
            df_linked_peaks_AT_temp = pd.merge_asof( df_linked_peaks_nucleotide, self.df_genome_seqs_duplicated_5X_grouped_AT_GTF_within, 

                        left_index=False, right_index=False, 
                        left_on="peak_location",  right_on="first", 
                        by=["seqname", "strand" ],
                        suffixes=('','_AT'), tolerance=tolerance, 
                        allow_exact_matches=True, direction=method)
            df_linked_peaks_AT = pd.concat([df_linked_peaks_AT, df_linked_peaks_AT_temp])

        df_linked_peaks_AT.drop_duplicates(subset =[ "seqname", "peak_location", "first", "file"], inplace=True)
        df_linked_peaks_AT.dropna(how ="any", subset =["seqname","peak_location", "first"], inplace=True)
        
        df_linked_peaks_AT.rename(
            columns = dict(zip([col for col in df_linked_peaks_AT.columns\
                                if col not in self.df_peak_links.columns and "AT" not in col], 
                                [col+"_AT" for col in df_linked_peaks_AT.columns\
                                if col not in self.df_peak_links.columns and "AT" not in col])),
            inplace=True)

        return df_linked_peaks_AT

    def AT_Stretch(self, df_linked_peaks_AT):
        df_peak_links = self.df_peak_links



        # What to do when the peak can be explained by an A/T stretch of a given length
        # assumption: longer A/T stretches are more likely to be primed --> tolrance for these peaks is therefore higher
        df_linked_peaks_AT["peak_to_geneUP_AT"] = df_linked_peaks_AT.apply(lambda row: abs(row["start_AT"] - row["peak_location"])if row["strand"] == "+" else abs(row["end_AT"] - row["peak_location"]), axis = 1)

        # identify AT-stretches
        pos_AT_unlinked = (df_linked_peaks_AT["strand"] == "+") &(df_linked_peaks_AT["start_preceding_gene_AT"] == df_linked_peaks_AT["start"]) & (df_linked_peaks_AT["start_AT"] != df_linked_peaks_AT["start"])
        neg_AT_unlinked = (df_linked_peaks_AT["strand"] == "-") &(df_linked_peaks_AT["end_next_gene_AT"] == df_linked_peaks_AT["end"]) & (df_linked_peaks_AT["end_AT"] != df_linked_peaks_AT["end"])

        # count_AT is the length of the A/T stretch
        count_AT =df_linked_peaks_AT["count_AT"]
        
        # AT_gene_buffer: distance from start of gene, AT_gene_buffer_factor: factor relative to the length of A/T stretch
        AT_gene_buffer, AT_gene_buffer_factor = 600, 4 
        AT_UPlimit = np.round(AT_gene_buffer/(np.floor(count_AT/AT_gene_buffer_factor)))
        peak_to_geneUP_AT_unlinked = (df_linked_peaks_AT["peak_to_geneUP_AT"] <= AT_UPlimit)

        # which peaks should not be considered because of their A/T stretch association
        df_unlinked_peaks_AT = df_linked_peaks_AT[(pos_AT_unlinked|neg_AT_unlinked) 
                                                & (peak_to_geneUP_AT_unlinked)
                                                ]["index"]

        # which peaks are identified as having AT-stretch association, and append info the the peaks dataframe
        AT_stretch = (df_peak_links["index"].isin(list(df_unlinked_peaks_AT)))
        self.df_peak_links["AT_Stretch"] = AT_stretch


        

    def blacklist(self, df_linked_peaks_AT):
        df_peak_links = self.df_peak_links
        df_linked_peaks_AT = pd.merge(df_summary, df_linked_peaks_AT, on = "file")

        # blacklist genes that have highscoring peaks in range of AT_rich region
        # these genes should be excluded from determining peak boundaries
        
        # identify peaks downstream of AT-stretches that fall within the next tandem gene
        black_list = ((df_linked_peaks_AT["score_peak"] > df_linked_peaks_AT["q90"]) 
                    & (((df_linked_peaks_AT["start"] == df_linked_peaks_AT["start_preceding_gene_AT"]) & (df_linked_peaks_AT["strand"] == "+") )
                    | ((df_linked_peaks_AT["start"] == df_linked_peaks_AT["start_next_gene_AT"]) & (df_linked_peaks_AT["strand"] == "-") ))
                    
                    & (df_linked_peaks_AT.duplicated(subset = ["file", "seqname", "strand", "start"], keep =False)))
        

        # create a dataframe with linking blacklisted peaks to AT-stretch dataframe 
        df_linked_peaks_AT_blacklist = pd.DataFrame()
        for strand, next_gene_start in zip(["+", "-"], ["start_preceding_gene_AT", "start_next_gene_AT"]):
            df_linked_peaks_AT_blacklist_temp = df_linked_peaks_AT[(black_list) & (df_linked_peaks_AT["strand"] == strand)][["file", "seqname", "strand",
                                                                                                                            "start_AT", "end_AT",
                                                                                                                            "start_preceding_gene_AT", "end_preceding_gene_AT", 
                                                                                                                            "start_next_gene_AT", "end_next_gene_AT",
                                                                                                                            "peak_location"]]
            df_linked_peaks_AT_blacklist_temp["start"] = df_linked_peaks_AT_blacklist_temp[next_gene_start]
            df_linked_peaks_AT_blacklist_temp["blacklist"] = True
            df_linked_peaks_AT_blacklist = pd.concat([df_linked_peaks_AT_blacklist_temp, df_linked_peaks_AT_blacklist])


        # merge the blacklist-AT dataframe with the comeplete peaks dataframe
        df_peak_links["strand_sign"] = df_peak_links.apply(lambda row: 1 if row["strand"]== "+" else -1, axis =1)
        df_peak_links_blacklist = pd.merge(df_peak_links, df_linked_peaks_AT_blacklist,  on = ["file", "seqname", "strand", "start"], suffixes=("", "_blacklist"), how = "outer")
        df_peak_links_blacklist.fillna({"blacklist":False, "peak_location_blacklist":df_peak_links_blacklist["peak_location"] - df_peak_links_blacklist["strand_sign"]*100 }, inplace=True)

        # identify blacklisted peaks that fall downstream of AT-stretch associated peak
        df_peak_links_blacklist_pos = df_peak_links_blacklist[
            ((df_peak_links_blacklist["blacklist"]) & ((df_peak_links_blacklist["strand"] == "+") & (df_peak_links_blacklist["peak_location"] >=  df_peak_links_blacklist["peak_location_blacklist"])))
            ]
        df_peak_links_blacklist_neg = df_peak_links_blacklist[
            ((df_peak_links_blacklist["blacklist"]) & ((df_peak_links_blacklist["strand"] == "-") & (df_peak_links_blacklist["peak_location"] <= df_peak_links_blacklist["peak_location_blacklist"])))
            ]
        df_peak_links_blacklist = pd.concat([df_peak_links_blacklist_pos, df_peak_links_blacklist_neg])

        # merge the blacklist dataframe with the comeplete peaks dataframe
        self.df_peak_links = pd.merge(df_peak_links, df_peak_links_blacklist[["file", "seqname", "strand", "start","blacklist",  "peak_location", "peak_location_blacklist"]],
                on = ["seqname", "strand", "start", "peak_location"], suffixes=("", "_bl"), how = "left").drop_duplicates(["file", "seqname", "strand", "start","blacklist",  "peak_location"], keep = "first")
        self.df_peak_links.fillna({"blacklist":False, "peak_location_blacklist":df_peak_links["peak_location"] - df_peak_links["strand_sign"]*100 }, inplace=True)

        
    
    def Clean(self):
        
        df_linked_peaks_AT = self.Link_w_AT()
        self.AT_Stretch(df_linked_peaks_AT)
        self.blacklist(df_linked_peaks_AT)

        return self.df_peak_links

class CollapsePeaks(PeakStreamClean):
    
    def __init__(self, bed_directory, coordinate_file, 
                 genome_fasta_dir, n, stretch_length,
                 new_coordinates_file_directory=None, coordinates_file_prefix = None, coordinates_file_suffix = None):
        
        PeakStreamClean.__init__(self,bed_directory, coordinate_file,
                                 genome_fasta_dir, n, stretch_length)

        self.df_peak_links = self.Clean()
        self.coordinate_file = coordinate_file
        self.bed_directory = bed_directory
        self.new_coordinates_file_directory = new_coordinates_file_directory
        self.coordinates_file_prefix = coordinates_file_prefix
        self.coordinates_file_suffix = coordinates_file_suffix
        

    def putativeUTRS(self):

        # function to define "UTRs" longer than 500bp
        for gene_position in ["start", "end"]:
            self.df_peak_links[f"peak_to_{gene_position}"] = self.df_peak_links["peak_location"] - self.df_peak_links[gene_position]
        self.df_peak_links["UTRgreater500"] = (((self.df_peak_links["strand"]== "-") & (self.df_peak_links["peak_to_start"] < -500)) 
                                               |
                                                 ((self.df_peak_links["strand"]== "+") & (self.df_peak_links["peak_to_end"] > 500)))
        
    def AT_stretch(self):

        # include peaks not in AT-strech - if peak is linked AT stretch and the it must not be part of long UTR
        self.df_peak_links = self.df_peak_links[(~self.df_peak_links["AT_Stretch"]) | ~((self.df_peak_links["AT_Stretch"] & (self.df_peak_links["UTRgreater500"])))]
    
    def countpeaks(self):
        # count number of peaks linked to single gene across all files
        self.df_peak_links["peak_count_per_gene"] = self.df_peak_links.groupby(["seqname","strand","start", "end", "peak_location"])["peak_location"].transform('count')

        # count number of peaks linked to single gene within each file
        self.df_peak_links["peaks_per_gene_per_file"] = self.df_peak_links.groupby(["file","seqname","strand", "start", "end"])["start"].transform('count')
    
    def aggregateLinks(self):

        # function to aggregate statistics on scores for each gene
        df_peak_links_groupby = self.df_peak_links.groupby(["file", "seqname","strand","start", "end"]).agg({"score_peak":["max", "median"]}).reset_index().T.reset_index()
        df_peak_links_groupby["column"] = df_peak_links_groupby["level_0"] + df_peak_links_groupby["level_1"]
        df_peak_links_groupby = df_peak_links_groupby.set_index("column").T.drop(["level_0", "level_1"])
        self.df_peak_links = pd.merge(self.df_peak_links, df_peak_links_groupby, on = ["file", "seqname","strand","start", "end"], how = "outer")

    def summary(self):
        # df_summary["q75"] to serve as internal control for noise in data
        self.df_peak_links = pd.merge(df_summary, self.df_peak_links, on = "file")

    def Rules(self):

        # categorizing rows with True/False statements e.g. whether a peak is linked to one or more than one gene within a file
        only_1PeakperGene_per_file = (self.df_peak_links["peaks_per_gene_per_file"] == 1)
        more_1PeakperGene_per_file = (self.df_peak_links["peaks_per_gene_per_file"]> 1)
        score_morethan_medianscore_per_file = (((self.df_peak_links["score_peak"] >= self.df_peak_links["score_peakmedian"]) 
                                                & (self.df_peak_links["score_peakmedian"] > self.df_peak_links["q75"])) 
                                                | (self.df_peak_links["score_peak"] == self.df_peak_links["score_peakmax"]))
        return only_1PeakperGene_per_file, more_1PeakperGene_per_file, score_morethan_medianscore_per_file
    
    def applyRules(self):

        only_1PeakperGene_per_file, more_1PeakperGene_per_file, score_morethan_medianscore_per_file = self.Rules()

        # application of rules
        self.df_peak_links["peak_links_thresh"] = only_1PeakperGene_per_file | ((more_1PeakperGene_per_file) & (score_morethan_medianscore_per_file))
        df_peak_links_thresh= self.df_peak_links[self.df_peak_links["peak_links_thresh"]]

        return df_peak_links_thresh
    
    def aggregateThresh(self):
 
        df_peak_links_thresh = self.applyRules()
        # function to aggregate statistics on distnace of peak to the end of the gene 
            # (negative stranded genes end at "start", positive end at "end")
        df_peak_links_thresh_groupby = df_peak_links_thresh.groupby(["file", "seqname","strand","start", "end"]).agg({"peak_to_start":["max", "min"], "peak_to_end":["max", "min"]}).reset_index().T.reset_index()
        df_peak_links_thresh_groupby["column"] = df_peak_links_thresh_groupby["level_0"] + df_peak_links_thresh_groupby["level_1"]
        df_peak_links_thresh_groupby = df_peak_links_thresh_groupby.set_index("column").T.drop(["level_0", "level_1"])
        self.df_peak_links = pd.merge(self.df_peak_links, df_peak_links_thresh_groupby, on = ["file", "seqname","strand","start", "end"], how = "outer")

    def Collapse(self):

        # putting everything together to define gene boundaries based on transcript data
        self.df_peak_links["last_peak_link"] = ((self.df_peak_links["strand"] == "+") & (self.df_peak_links["peak_to_end"] == self.df_peak_links["peak_to_endmax"]) 
                                                |
                                                  (self.df_peak_links["strand"] == "-") & (self.df_peak_links["peak_to_start"] == self.df_peak_links["peak_to_startmin"]))
        self.df_peak_links["UTRgreater500"] = (((self.df_peak_links["strand"]== "-") & (self.df_peak_links["peak_to_start"] < -500)) 
                                               |
                                                 ((self.df_peak_links["strand"]== "+") & (self.df_peak_links["peak_to_end"] > 500)))
        self.df_peak_links["possible_novel_transcript"] = (self.df_peak_links["UTRgreater500"] ) & (~(self.df_peak_links["last_peak_link"]) | (self.df_peak_links["blacklist"]))
        self.df_peak_links["novel_transcript"] = (self.df_peak_links["UTRgreater500"]) & ((self.df_peak_links["last_peak_link"]) | (self.df_peak_links["blacklist"]))
        false_pos_neg = (self.df_peak_links["strand"] == "-") & (self.df_peak_links["peak_to_start"]<= -500) & (self.df_peak_links["novel_transcript"]==False) & (self.df_peak_links["peak_to_end"]<0)
        false_pos_pos = (self.df_peak_links["strand"] == "+") & (self.df_peak_links["peak_to_end"]>= 500) & (self.df_peak_links["novel_transcript"]==False) & (self.df_peak_links["peak_to_start"]>0)
        self.df_peak_links["False_pos"] = false_pos_neg|false_pos_pos
        self.df_peak_links = self.df_peak_links[self.df_peak_links["False_pos"] == False]
        self.df_peak_links["peak_count_per_transcript"] = self.df_peak_links.groupby(["seqname","strand","start", "end", "possible_novel_transcript", "peak_location" ])["peak_location"].transform('count')
        self.df_peak_links["transcript_peak"]  = (self.df_peak_links["peak_links_thresh"] | self.df_peak_links["novel_transcript"] )& (self.df_peak_links["peak_count_per_transcript"] > 2)

        self.df_peak_links = self.df_peak_links[self.df_peak_links["transcript_peak"] ]

    def MatchtoGTF(self):
        # match peaks back to original GTF file

        df_peak_links_GTF_GeneRef = pd.merge(self.dfGTF_GeneRef, self.df_peak_links, on = ['seqname','strand','start','end'], how = "outer", suffixes=("", "peaks"))

        df_peak_links_GTF_GeneRef['peak_location'].fillna(df_peak_links_GTF_GeneRef['end'], inplace=True)
        df_peak_links_GTF_GeneRef['score_peak'].fillna(-20, inplace=True)
        df_peak_links_GTF_GeneRef['novel_transcript'].fillna(False, inplace=True)
        df_peak_links_GTF_GeneRef['blacklist'].fillna(False, inplace=True)

        df_peak_links_GTF_GeneRef["no_peak"] =( df_peak_links_GTF_GeneRef["score_peak"] == -20)
        df_peak_links_GTF_GeneRef[["peak_to_start", "peak_to_end"]].fillna(0, inplace=True)


        df_peak_links_GTF_GeneRef.drop_duplicates(["seqname", "strand","start", "end"], keep = "first")

        return df_peak_links_GTF_GeneRef
    
    def renameAttribute(self):
        # attribute renaming to include new transcript attribute
        df_peak_links_GTF_GeneRef = self.MatchtoGTF()

        df_peak_links_GTF_GeneRef["old gene_id"] = df_peak_links_GTF_GeneRef.apply(lambda row: '''gene_id "'''  + row["new A"] +  '''"''' if "ntar" not in row["attribute"] else "Parent="+ row["new A"] , axis=1)

        update_col = "updated_gene_id"
        df_peak_links_GTF_GeneRef[update_col] = df_peak_links_GTF_GeneRef.apply(lambda row: "h" + row["new A"] if row["novel_transcript"] and not row["no_peak"] else row["new A"], axis = 1)
        df_peak_links_GTF_GeneRef[update_col] = df_peak_links_GTF_GeneRef.apply(lambda row: "bl" + row[update_col] if row["blacklist"] and row["novel_transcript"] and not row["no_peak"] else row[update_col], axis = 1)

        df_peak_links_GTF_GeneRef["new gene_id"] = df_peak_links_GTF_GeneRef.apply(lambda row: '''gene_id "'''  + row[update_col] +  '''"''' if row["new A"] in row["attribute"] else ''' gene_id "'''  + row[update_col] +  '''"''', axis=1)
        df_peak_links_GTF_GeneRef['new_attribute'] = df_peak_links_GTF_GeneRef.apply(lambda row: row["attribute"].replace(row['old gene_id'],row["new gene_id"]), axis=1)


        df_peak_links_GTF_GeneRef['new_attribute'] = df_peak_links_GTF_GeneRef.apply(lambda row: "ID="+row[update_col]+";Name="+row[update_col], axis=1)

        return df_peak_links_GTF_GeneRef

    def value(self, x):
        return np.unique(x)[0]

    
    
    def aggregatorGTF(self):
        df_peak_links_GTF_GeneRef = self.renameAttribute()
        aggragator = {"peak_location":["min", "max"], "score_peak": ["mean"], "attribute": self.value }
        df_new_GTF = df_peak_links_GTF_GeneRef.groupby(["seqname", "source", "feature","strand", "start", "end", "score", "frame", "new_attribute", 
                                                        "novel_transcript", "blacklist", "no_peak",]).agg(aggragator).reset_index().T.reset_index()
        df_new_GTF["column"] = df_new_GTF["level_0"] + df_new_GTF["level_1"]
        df_new_GTF = df_new_GTF.set_index("column").T.drop(["level_0", "level_1"])

        for gene_position, m in zip(["start", "end"], ["min", "max"]):
            df_new_GTF[f"peak_to_{gene_position}"] = df_new_GTF[f"peak_location{m}"] - df_new_GTF[gene_position]

        
        return df_new_GTF

    def hypothetical_correction(self, df_new_GTF):
        # some genes were overwritten as hypothetical, this code adds a new line for each overwritten gene
        hypothetical_only = df_new_GTF[(df_new_GTF["new_attribute"].str.contains("h")) & (~df_new_GTF["attributevalue"].duplicated())]
        hypothetical_only[["novel_transcript","blacklist","no_peak", "score_peak"]] = False,False, True, -20
        
        hypothetical_only["peak_locationmin"] = hypothetical_only["start"]
        hypothetical_only["peak_locationmax"]= hypothetical_only["end"]
        hypothetical_only["new_attribute"] = hypothetical_only["new_attribute"].str.replace("h", "").str.replace("bl", "")
        df_new_GTF = df_new_GTF.append(hypothetical_only)

        return df_new_GTF

    
    def plot_result(self, df_new_GTF):

        fig, axes =plt.subplots(figsize = (10,5), ncols = 2, sharey= True)
        markers = {True: "X", False: "o"}
        dict_strand = dict(zip(["-", "+"], [-1, 1]))
        
        for ax, peak_to, strand in zip(axes, ["start", "end"], ["-", "+"]):
            data = df_new_GTF[df_new_GTF["strand"] ==strand]
            data["blacklist_actual"] = (data["blacklist"]) & (abs(data[f"peak_to_{peak_to}"]) > 500)


            sns.scatterplot(data = data, x =f"peak_to_{peak_to}", y = "score_peakmean", ax =ax, style = "novel_transcript",
                            hue = "blacklist_actual", hue_order = [False, True],  palette = ["#012536", "#f04175"],markers =markers,
                            s= 20, legend = True, alpha = .5)
            ax.set_yscale("log")
            for x in [1000,500,-400]:
                x= dict_strand[strand]*x
                ax.vlines(ymin = 1, ymax = 10**7, x = x, color = "k", lw = 0.5)
                ax.text(x = x, y = 10**7, s = f"{x}", fontsize=6)
            ax.hlines(xmin = ax.get_xlim()[0], xmax = ax.get_xlim()[1], y = np.mean(df_summary["q75"]), color = "k", lw = 0.5)
            dict_strand = dict(zip(["-", "+"], [-1, +1]))

            ax.set_xlim(np.sort((dict_strand[strand]*5000, dict_strand[strand]*-1000)))

        plt.show()

    def new_start(self, row):
        updated_end = min(600, abs(row["start"]-row["end"]))
        if row["novel_transcript"] == True:
            new_start = row["peak_locationmax"] - 100
            return max(1, new_start)

        elif row["strand"] == "-":
            new_start = min(row["start"], row["peak_locationmin"]) - 100
            return max(1, new_start)
        else:
            new_start = row["end"] - updated_end
            return max(1, new_start)

    def new_end(self, row):
        updated_end = min(600, abs(row["start"]-row["end"]))
        if row["novel_transcript"] == True:
            new_end = row["peak_locationmax"] + 100
            return max(1, new_end)

        elif row["strand"] == "+":
            new_end =max(row["end"], row["peak_locationmax"]) + 100
            return max(1, new_end)
        else:
            new_end = row["start"]+ updated_end
            return max(1, new_end)
        
    def CollapseGTF(self):
        df_new_GTF = self.aggregatorGTF()
        df_new_GTF = self.hypothetical_correction(df_new_GTF)
        df_new_GTF.sort_values(["seqname", "start", "end"], inplace=True)
        df_new_GTF.rename(columns  = {"new_attribute": "attribute"}, inplace=True)
        df_new_GTF = df_new_GTF[df_new_GTF["feature"] == "CDS"]
        df_new_GTF["new_start"]  = df_new_GTF.apply(lambda row: self.new_start(row) , axis =1)
        df_new_GTF["new_end"]  = df_new_GTF.apply(lambda row: self.new_end(row) , axis =1)
        df_new_GTF["start_diff"] = df_new_GTF["new_start"] -  df_new_GTF["start"]
        df_new_GTF["end_diff"] = df_new_GTF["new_end"] -  df_new_GTF["end"]
        df_new_GTF["gene_size"] = abs(df_new_GTF["start"] -  df_new_GTF["end"])
        
        return df_new_GTF

    def MakeSaveGTF(self):
        df_new_GTF = self.CollapseGTF()
        df_new_GTF["start"] = df_new_GTF["new_start"].astype(int)
        df_new_GTF["end"] = df_new_GTF["new_end"].astype(int)
        df_new_GTF["frame"] = "."
        df_new_GTF["feature"] = "predicted_3primeUTR"
        df_new_GTF["source"] = "PeakStream"

        df_new_GTF_save = df_new_GTF[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]
    
        return df_new_GTF_save
    
    def outdir(self):
        if not self.new_coordinates_file_directory:
            # infer the directory name from the input folder
            self.new_coordinates_file_directory = os.path.join(self.bed_directory.split("bedgraph")[0], "PeakStreamOutput")

        if not os.access(self.new_coordinates_file_directory, os.W_OK):
            # if no write persmissions, write file to current directory in PeakStreamOutput directory
            self.new_coordinates_file_directory = "PeakStreamOutput"
        
        if not os.path.exists(self.new_coordinates_file_directory):
            os.mkdir(self.new_coordinates_file_directory)

    
    def prefix(self):
        if not self.coordinates_file_prefix:
            self.coordinates_file_prefix = self.coordinate_file.split("/")[-1].split(".")[0]

    def suffix(self):
        if not self.coordinates_file_suffix:
            self.coordinates_file_suffix = [i for i in self.bed_directory.split("bedgraph")[0].split("/") if i][-1]


    
    def save_pygenometracks(self, df_new_GTF_save, new_coordinates_file_path):

        pygenome_coordinates_file_path = new_coordinates_file_path.split(".")[0]+"_pygenometracks.gff"
        df_new_GTF_save["feature"] = "exon"
        df_GTF = self.dfGTF_GeneRef[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand','frame', 'attribute']]

        df_GTF["end"] = df_GTF.apply(lambda row: row["end"]-min(600, abs(row["start"] - row["end"])) if row["strand"] == "+" and row["feature"] == "CDS" else row["end"], axis=1)
        df_GTF["start"] = df_GTF.apply(lambda row:  row["start"]+ min(600, abs(row["start"] - row["end"])) if row["strand"] == "+" and row["feature"] == "CDS" else row["start"], axis=1)
        df_new_GTF_pygenome = pd.concat([df_new_GTF_save, df_GTF])
        df_new_GTF_pygenome.to_csv(pygenome_coordinates_file_path, sep = "\t", index=False, header = None,  quoting=csv.QUOTE_NONE)




    
    def SaveNewGTF(self):
        self.putativeUTRS()
        self.AT_stretch()
        self.countpeaks()
        self.aggregateLinks() 
        self.summary()
        self.Rules()
        self.applyRules()
        self.aggregateThresh()
        self.Collapse()

        # # ## change this back later
        # df_new_GTF = self.aggregatorGTF()
        # return df_new_GTF


        df_new_GTF_save = self.MakeSaveGTF()
        self.prefix() #infer prefix bed-directory name if none was given
        self.suffix() #infer suffix bed-directory name if none was given
        self.outdir() # check out directory, make a new one if none was specified
        
        new_coordinates_file_path = os.path.join(self.new_coordinates_file_directory, f"{self.coordinates_file_prefix}_{self.coordinates_file_suffix}.gff")
        df_new_GTF_save.to_csv(new_coordinates_file_path, sep = "\t", index=False, header = None,  quoting=csv.QUOTE_NONE)
        self. save_pygenometracks(df_new_GTF_save, new_coordinates_file_path)

        return new_coordinates_file_path



def peakstream(bed_directory, coordinate_file, 
                 genome_fasta_dir, n, stretch_length,
                 new_coordinates_file_directory=None,
                 coordinates_file_prefix = None,
                 coordinates_file_suffix = None, 
                 ):
    
    """
    Executes the PeakStream pipeline to process bedgraph files and 
    identify transcript boundaries based on poly-A priming RNA-seq data.

    Parameters:
    - bed_directory (str): Path to the directory containing bedgraph files 
                           (signal intensity data).
    - coordinate_file (str): Path to the GTF file containing gene annotations.
    - genome_fasta_dir (str): Path to the directory containing genome FASTA files.
    - n (int): Number of data points to smooth or window size for peak detection.
    - stretch_length (int): Minimum length of A/T stretches to identify.
    - new_coordinates_file_directory (str, optional): Directory to save the new 
                                                      GTF file with updated 
                                                      transcript coordinates. 
                                                      Default is None.
    - coordinates_file_prefix (str, optional): Prefix for the new GTF file name.
                                               Default is None.
    - coordinates_file_suffix (str, optional): Suffix for the new GTF file name.
                                               Default is None.

    Returns:
    - str: Path to the newly saved GTF file containing refined transcript coordinates.
    """


    coordinates_path = CollapsePeaks(bed_directory, coordinate_file,  genome_fasta_dir, n, stretch_length,
                 new_coordinates_file_directory, coordinates_file_prefix, coordinates_file_suffix, 
                 ).SaveNewGTF()
    
    return coordinates_path