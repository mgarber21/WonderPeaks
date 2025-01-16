import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Parent_directory = os.path.dirname(os.getcwd())  # Parent file's directory
# sys.path.insert(0, Parent_directory)
from .Helper_Functions import *


class ATstretches():
    def __init__(self, directory, coordinate_file, genome_fasta_dir, stretch_length):
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
        self.gtf_data  = load_coordinates(coordinate_file = coordinate_file)
        self.directory = directory
        self.PS_directory = PeakStream_Direcory(self.directory)


    def find_AT_stretches(self):
        """
        Identify stretches of A/T nucleotides exceeding a defined length.

        Returns:
        - pd.DataFrame: DataFrame of A/T stretches exceeding the specified length.
        """
        print("Searching for A/T stretches within genome...")

        stretch_coordinate_file = os.path.join(
            self.PS_directory,
            os.path.basename(self.genome_fasta_dir).replace(".fasta", "_ATstretch_coordinates.csv")
            )

        if os.path.isfile(stretch_coordinate_file):
            return pd.read_csv(stretch_coordinate_file, index_col = False)
        
        with open(self.genome_fasta_dir) as fasta_file:  # Will close handle cleanly
            identifiers = []
            seqs = []
            iter = []
            for title, sequence in SimpleFastaParser(fasta_file):
                identifiers.append(title.split(None, 1)[0])  # First word is ID
                seqs.append(list(sequence))
                iter.append(list(range(len(sequence))))

        
        # iterate over nucleotide sequences and make a dataframe of runs of nucleotides
        genome_df = pd.DataFrame([identifiers, seqs, iter]).T.explode([1, 2]).reset_index()
        genome_df.rename(columns={0:"seqname", 1:"nucleotide", 2:"iter"}, inplace=True)
        genome_df['nucleotide_run'] = (genome_df['nucleotide'] != genome_df['nucleotide'].shift(1)).cumsum()
        genome_2X_df = genome_df[genome_df["nucleotide_run"].duplicated(keep = False)]
        subgroup_counts = genome_2X_df['nucleotide_run'].value_counts()

        # 
        filtered = genome_2X_df[genome_2X_df['nucleotide_run'].isin(subgroup_counts[subgroup_counts>=self.stretch_length].index)]
        filtered_grouped = filtered.groupby(["seqname", "nucleotide", "nucleotide_run"],as_index=False).agg({"iter":["first", "last","count"]}).droplevel(axis=1, level = 0)
        filtered_grouped.columns = ["seqname", "nucleotide", "nucleotide_run"] +['first', 'last', 'count']
        filtered_grouped_AT= filtered_grouped[filtered_grouped["nucleotide"].isin(["A", "T"])]

        filtered_grouped_AT.to_csv(stretch_coordinate_file, index = False)
    
        return filtered_grouped_AT
    
    def plot_AT_dist(self):
        """
        Plot distribution of stretches of A/T nucleotides exceeding a defined length.

        Returns:
        - a plot of the distribution of A/T stretches with your genome.
        """
        filtered_grouped_AT = self.find_AT_stretches()
        fig, ax = plt.subplots(figsize = (15,5))
        sns.barplot(data = pd.DataFrame(filtered_grouped_AT.value_counts(["count", "nucleotide"])).reset_index(),
                    x="count",y=0,
                    hue = "nucleotide", ax = ax, palette= ["#000000", "#000587"],
                    linewidth=2.5, edgecolor="#ffffff"
                    )
        ax.set_yscale("log")
        ax.set_ylabel("# of A of T stretches of {x} length")
        ax.set_xlabel("Length of A or T stretches")

        ax.bar_label(ax.containers[0], fontsize=10, padding = 1, color ="#000000" )
        ax.bar_label(ax.containers[1], fontsize=10, padding = 10, color ="#000587")
        ax.set_xlim(-1, 27)
        ax.set_ylim(ax.get_ylim()[0], 10**5)
        ax.set_title("Distribution of A or T stretches across the genome")
        ax.legend(loc = "upper right", ncol =2, frameon =False)
        sns.despine(top =True, right = True)
        plt.show()
    
    def merge_AT_w_coordinates(self):
        """
        Link stretches of A/T nucleotides exceeding a defined length to the coordinates of genes
        within your genome

        Returns:
        - pd.DataFrame: DataFrame of A/T stretches exceeding the specified length merged with your genomes coordinates
        """

        print("Linking A/T stretches within genome to gene coordinates...")
        filtered_grouped_AT = self.find_AT_stretches()


        filtered_grouped_AT['first'] = filtered_grouped_AT['first'].astype("int")
        self.gtf_data['start'] = self.gtf_data['start'].astype("int")
        filtered_grouped_AT.sort_values(by ="first", inplace=True)
        self.gtf_data.sort_values(by ="start", inplace=True)
        filtered_grouped_AT["strand"] = filtered_grouped_AT.apply(lambda row: "+" if row["nucleotide"] == "A" else "-", axis= 1)


        AT_coordinates = pd.DataFrame()


        for method in ["backward", "forward"]:
            self.gtf_data["start"] = self.gtf_data["start"].astype("int")
            self.gtf_data.sort_values(by ="start", inplace=True)
            AT_coordinates_temp = pd.merge_asof(filtered_grouped_AT, self.gtf_data[self.cols], 
                        left_index=False, right_index=False, 
                        left_on=["first"], right_on="start", 
                        by=["seqname", "strand"],
                        suffixes=('_x', '_y'), tolerance=20000, 
                        allow_exact_matches=True, direction=method)
            AT_coordinates = pd.concat([AT_coordinates, AT_coordinates_temp])

        AT_coordinates.dropna(how ="any", inplace=True)
        AT_coordinates.drop_duplicates(subset =["seqname", "first", "last", "start", "end"], inplace=True)

        return AT_coordinates
    

    def ApplyATStretches(self):

        """
        Identify stretches of A/T nucleotides exceeding a defined length that lie within a gene in you genome

        Returns:
        - pd.DataFrame: DataFrame of A/T stretches exceeding the specified length linked to genes in your genome
        """

        AT_coordinates = self.merge_AT_w_coordinates()
        # AT regions witin the first 85% of a gene if the gene >

        #gene:             ------>
        #ATcoverage:       -----
        pos_AT = ((AT_coordinates["strand"] == "+"))
        neg_AT = ((AT_coordinates["strand"] == "-"))

        AT_coordinates_inGene =AT_coordinates[pos_AT|neg_AT]
        inGene_gene = (
            (AT_coordinates_inGene["first"] > AT_coordinates_inGene["start"])
                & (AT_coordinates_inGene["first"] < AT_coordinates_inGene["end"])
                )
        AT_coordinates_inGene = AT_coordinates_inGene[inGene_gene]

        return AT_coordinates_inGene
    

