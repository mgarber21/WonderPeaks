# updates after new year
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


from PeakStream.Helper_Functions import *
from PeakStream.AT_stretches import *
from PeakStream.WonderPeaks4UTRs import *
from PeakStream.PeakLinks import *
from PeakStream.PeakFlags import *
from PeakStream.feature_counts import *



class PeakStream(PeakStreamFlag):
    """
    Final processing of PeakStream    
    
    """
    
    def __init__(self, directory, coordinate_file, 
                 genome_fasta_dir, n, stretch_length,
                 User_inputs_dict,
                 new_coordinates_file_directory=None, 
                 coordinates_file_prefix = None, coordinates_file_suffix = None,
                 gene_biotypes = "protein_coding",
                 rerun = False
                 ):
        
        
        
        
        self.directory = directory
        self.PS_directory = os.path.join(self.directory, "PeakStream")
        self.coordinate_file = coordinate_file
        self.new_coordinates_file_directory = new_coordinates_file_directory
        self.coordinates_file_prefix = coordinates_file_prefix
        self.coordinates_file_suffix = coordinates_file_suffix
        self.prefix = prefix(self.coordinate_file, self.coordinates_file_prefix) #infer prefix directory name if none was given
        self.suffix = suffix(self.PS_directory, self.coordinates_file_suffix) #infer suffix directory name if none was given
        self.outdir = outdir(self.PS_directory, self.new_coordinates_file_directory) # check out directory, make a new one if none was specified
        self.new_coordinates_file_path = os.path.join(
                                                self.outdir, 
                                                f"{self.prefix}_{self.suffix}.gtf"
                                                )
        self.PS_data_output = os.path.join(self.PS_directory, "PeakStream_new_coordinates_data.csv")
        self.User_inputs_dict = User_inputs_dict
        self.designfactor = self.User_inputs_dict["designfactor"]
        self.rerun = rerun
        self.gene_biotypes = [i.strip(" ") for i in gene_biotypes.split(";")]


        # Check if the output file already exists
        if not os.path.isfile(self.new_coordinates_file_path) or self.rerun:
            PeakStreamFlag.__init__(self,directory, coordinate_file, 
                genome_fasta_dir,  n, stretch_length)
            self.peak_links = self.ApplyFlags()
        else:
            self.df_summary = pd.read_csv(os.path.join(self.PS_directory, "bedgraph_summary.csv"))
            self.new_coordinates_GTF = pd.read_csv(self.new_coordinates_file_path, sep = "\t", header = None)
            self.new_coordinates = pd.read_csv(self.PS_data_output)
            self.gtf_data  = load_coordinates(coordinate_file = coordinate_file)

            
        
        

    def putativeUTRS(self):
        """
        Includes rows if a peak is a possible UTR 

        """

        # function to define "UTRs" longer than 500bp
        for gene_position in ["start", "end"]:
            self.peak_links[f"peak_to_{gene_position}"] = self.peak_links["peak_location"] - self.peak_links[gene_position]
        self.peak_links["UTRgreater500"] = (((self.peak_links["strand"]== "-") & (self.peak_links["peak_to_start"] < -500)) 
                                               |
                                                 ((self.peak_links["strand"]== "+") & (self.peak_links["peak_to_end"] > 500)))
        
    def AT_stretch(self):
        """
        Includes rows if peak is not associated with AT stretch 

        """

        # include peaks not in AT-strech - if peak is linked AT stretch and the it must not be part of long UTR
        self.peak_links = self.peak_links[(~self.peak_links["AT_Stretch"]) | ~((self.peak_links["AT_Stretch"] & (self.peak_links["UTRgreater500"])))]
    
    def countpeaks(self):

        """
        counts number of peaks associated with a gene
        counts number of peaks associated with a gene and file

        """
        # count number of peaks linked to single gene across all files
        self.peak_links["peak_count_per_gene"] = self.peak_links.groupby(["seqname","strand","start", "end", "peak_location"])["peak_location"].transform('count')

        # count number of peaks linked to single gene within each file
        self.peak_links["peaks_per_gene_per_file"] = self.peak_links.groupby(["file","seqname","strand", "start", "end"])["start"].transform('count')
    
    def aggregateLinks(self):
        """
        aggregate statistics on scores for each gene

        """

        # function to aggregate statistics on scores for each gene
        peak_links_groupby = self.peak_links.groupby(["file", "seqname","strand","start", "end"]).agg({"score_peak":["max", "median"]}).reset_index().T.reset_index()
        peak_links_groupby["column"] = peak_links_groupby["level_0"] + peak_links_groupby["level_1"]
        peak_links_groupby = peak_links_groupby.set_index("column").T.drop(["level_0", "level_1"])
        self.peak_links = pd.merge(self.peak_links, peak_links_groupby, on = ["file", "seqname","strand","start", "end"], how = "outer")

    def applySummary(self):
        """
        apply internal control for noise in data

        """

        # self.df_summary["q75"] to serve as internal control for noise in data
        self.peak_links = pd.merge(self.df_summary, self.peak_links, on = "file")

    def setRules(self):
        """
        set rules on whether a peak is a UTR

        """

        # categorizing rows with True/False statements e.g. whether a peak is linked to one or more than one gene within a file
        only_1PeakperGene_per_file = (self.peak_links["peaks_per_gene_per_file"] == 1)
        more_1PeakperGene_per_file = (self.peak_links["peaks_per_gene_per_file"]> 1)
        score_morethan_medianscore_per_file = (((self.peak_links["score_peak"] >= self.peak_links["score_peakmedian"]) 
                                                & (self.peak_links["score_peakmedian"] > self.peak_links["q75"])) 
                                                | (self.peak_links["score_peak"] == self.peak_links["score_peakmax"]))
        return only_1PeakperGene_per_file, more_1PeakperGene_per_file, score_morethan_medianscore_per_file
    
    def applyRules(self):
        """
        apply set rules on whether a peak is a UTR

        """
        only_1PeakperGene_per_file, more_1PeakperGene_per_file, score_morethan_medianscore_per_file = self.setRules()

        # application of rules
        self.peak_links["peak_links_thresh"] = only_1PeakperGene_per_file | ((more_1PeakperGene_per_file) & (score_morethan_medianscore_per_file))
        peak_links_thresh= self.peak_links[self.peak_links["peak_links_thresh"]]

        return peak_links_thresh
    
    def aggregateThresh(self):
        """
        apply grouping to putative UTRs based on file

        """

        # This is where we have the trouble
        peak_links_thresh = self.applyRules()
        # function to aggregate statistics on distnace of peak to the end of the gene 
            # (negative stranded genes end at "start", positive end at "end")
        peak_links_thresh_groupby = peak_links_thresh.groupby(["file", "seqname","strand","start", "end"]).agg({"peak_to_start":["max", "min"], "peak_to_end":["max", "min"]}).reset_index().T.reset_index()
        peak_links_thresh_groupby["column"] = peak_links_thresh_groupby["level_0"] + peak_links_thresh_groupby["level_1"]
        peak_links_thresh_groupby = peak_links_thresh_groupby.set_index("column").T.drop(["level_0", "level_1"])
        self.peak_links = pd.merge(self.peak_links, peak_links_thresh_groupby, on = ["file", "seqname","strand","start", "end"], how = "outer")

    def Collapse(self):

        # putting everything together to define gene boundaries based on transcript data
        self.peak_links["last_peak_link"] = ((self.peak_links["strand"] == "+") & (self.peak_links["peak_to_end"] == self.peak_links["peak_to_endmax"]) 
                                                |
                                                  (self.peak_links["strand"] == "-") & (self.peak_links["peak_to_start"] == self.peak_links["peak_to_startmin"]))
        self.peak_links["UTRgreater500"] = (((self.peak_links["strand"]== "-") & (self.peak_links["peak_to_start"] < -500)) 
                                               |
                                                 ((self.peak_links["strand"]== "+") & (self.peak_links["peak_to_end"] > 500)))
        self.peak_links["possible_novel_transcript"] = (self.peak_links["UTRgreater500"] ) & (~(self.peak_links["last_peak_link"]) 
                                                                                              | (self.peak_links["blacklist"])
                                                                                              )
        # novel transcript is anything greater than 500bp away from gene
        self.peak_links["novel_transcript"] = (self.peak_links["UTRgreater500"]) | (self.peak_links["UTRgreater500"]) & ((self.peak_links["last_peak_link"]) 
                                                                                                                         | (self.peak_links["blacklist"])
                                                                                                                         )
        false_pos_neg = (self.peak_links["strand"] == "-") & (self.peak_links["peak_to_start"]<= -500) & (self.peak_links["novel_transcript"]==False) & (self.peak_links["peak_to_end"]<0)
        false_pos_pos = (self.peak_links["strand"] == "+") & (self.peak_links["peak_to_end"]>= 500) & (self.peak_links["novel_transcript"]==False) & (self.peak_links["peak_to_start"]>0)
        self.peak_links["False_pos"] = false_pos_neg|false_pos_pos
        self.peak_links = self.peak_links[self.peak_links["False_pos"] == False]
        self.peak_links["peak_count_per_transcript"] = self.peak_links.groupby(["seqname","strand","start", "end", "possible_novel_transcript", "peak_location" ])["peak_location"].transform('count')

        #  transcript peak is when last_peak_link is True
        self.peak_links["transcript_peak"]  = (self.peak_links["peak_links_thresh"] | self.peak_links["novel_transcript"] )& (self.peak_links["peak_count_per_transcript"] > 2 & (self.peak_links["last_peak_link"]))

        self.peak_links = self.peak_links[self.peak_links["transcript_peak"] ]

    def MatchtoGTF(self):
        # match peaks back to original GTF file

        peak_links2GTF = pd.merge(self.gtf_data, self.peak_links, on = ['seqname','strand','start','end'], how = "outer", suffixes=("", "peaks"))

        peak_links2GTF['peak_location'].fillna(peak_links2GTF['end'], inplace=True)
        
        # .fillna(peak_links2GTF['end'], inplace=True)
        peak_links2GTF['score_peak'].fillna(-20, inplace=True)
        peak_links2GTF['novel_transcript'].fillna(False, inplace=True)
        peak_links2GTF['blacklist'].fillna(False, inplace=True)

        peak_links2GTF["no_peak"] =( peak_links2GTF["score_peak"] == -20)
        peak_links2GTF[["peak_to_start", "peak_to_end"]].fillna(0, inplace=True)


        

        return peak_links2GTF
    
    def Attribute2GTF(self, peak_links2GTF, attribute_col = "attribute"):
        """
        attribute renaming to include new transcript attribute based on GTF format

        """ 

        # attribute renaming to include new transcript attribute based on GTF format:
        
        peak_links2GTF["gene_biotype"] = peak_links2GTF.apply(
                lambda row: rename_attributes(
                    row[attribute_col],  "gene_biotype")\
                        if rename_attributes(row[attribute_col],  "gene_biotype") != None\
                            else rename_attributes(row[attribute_col]),
                            axis=1)
        
        for named_attribute in ["gene_id", "gene_name"]:
            peak_links2GTF[named_attribute] = peak_links2GTF.apply(
                lambda row: rename_attributes(
                    row[attribute_col], named_attribute)\
                        if rename_attributes(row[attribute_col], named_attribute) != None\
                            else rename_attributes(row[attribute_col]),
                            axis=1)
            
            update_col = f"updated_{named_attribute}"
            peak_links2GTF[update_col] = peak_links2GTF.apply(
                lambda row: "h" + row[named_attribute].replace("h", "") if row["novel_transcript"] and not row["no_peak"] else row[named_attribute],
                axis = 1
                )
            # peak_links2GTF[update_col] = peak_links2GTF.apply(
            #     lambda row: "bl" + row[update_col].replace("bl", "") if row["blacklist"] and row["novel_transcript"] and not row["no_peak"] else row[update_col],
            #     axis = 1
            #     )


        # using GTF fomatting for attribute column
        peak_links2GTF['new_attribute'] = peak_links2GTF.apply(lambda row: f'''gene_id "{row["updated_gene_id"]}"; gene_name "{row["updated_gene_name"]}"; gene_biotype "{row["gene_biotype"]}";''',
                                                                            # "ID"+row["updated_gene_id"]+";Name="+row["updated_gene_name"], #GTF format
                                                                    axis=1)

        
        return peak_links2GTF

    def value(self, x):
        return np.unique(x)[0]

    
    
    def aggregatorGTF(self, peak_links2GTF):
        """
        aggregrate peaks based on attribute
        """
        
        aggragator = {"peak_location":["min", "max"], "score_peak": ["mean"], "attribute": self.value }
        new_coordinates = peak_links2GTF.groupby([
            "seqname", "source", "feature","strand", 
            "start", "end", "score", "frame", "new_attribute", 
            "novel_transcript", "no_peak", 
            "TandemGene", "dist_to_TandemGene", 
            # "blacklist"
            ]).agg(aggragator).reset_index().T.reset_index()
        new_coordinates["column"] = new_coordinates["level_0"] + new_coordinates["level_1"]
        new_coordinates = new_coordinates.set_index("column").T.drop(["level_0", "level_1"])

        for gene_position, m in zip(["start", "end"], ["min", "max"]):
            new_coordinates[f"peak_to_{gene_position}"] = new_coordinates[f"peak_location{m}"] - new_coordinates[gene_position]

        
        return new_coordinates

    
    def CDScounts_correction(self, new_coordinates):

        """
        take counts from featurecounts into account
        if counts across dataset for a given gene is less than 10, 
        and the peak is greater than 100bp away from the gene,
        CDScounts_correction will flag this gene as hypothetical

        """

        metadata =  metadata_upload(self.directory,self.designfactor,
                    meta_delim_whitespace = False,
                    meta_index_col = None)
        
        featurecounts_file = featurecounts(self.directory, self.User_inputs_dict, s=1, t="gene", T=20, 
                                        stardir="starout", featurecounts_subDir="featurecounts", rerun=self.rerun)
        featurecounts_file_biotype = featurecounts_biotype(featurecounts_file, self.gene_biotypes, t="gene")

        CDS_counts = pd.read_csv(
                    featurecounts_file_biotype,
                    )
        

        ID_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
        counts_cols = [col for col in CDS_counts.columns if col not in ID_cols]
        counts_basename = [get_basename(col) for col in counts_cols]

        CDS_counts.rename(
            columns=dict(zip(counts_cols, counts_basename)),
            inplace= True
            )

        CDS_counts_melt = pd.melt(CDS_counts, id_vars=ID_cols, var_name="basename",  value_vars=counts_basename, value_name = "count")
        CDS_counts_melt_merge = pd.merge(CDS_counts_melt, metadata, on = "basename")
        
        # gets the median count value across each designfactors
        CDS_counts_median = pd.DataFrame(CDS_counts_melt_merge.groupby(ID_cols+[self.designfactor])["count"].median()).reset_index()

        # gets the maximum median counts across all experiments for each gene
        CDS_counts_max_median = pd.DataFrame(CDS_counts_median.groupby(ID_cols)["count"].max()).reset_index()
        Geneids_count0 = CDS_counts_max_median[CDS_counts_max_median['count']< 10]["Geneid"]

        new_coordinates["Geneid"] = new_coordinates.apply(
            lambda row: rename_attributes(row["attributevalue"]), axis=1)

        new_coordinates["CDS_count0"] = (
            (new_coordinates["Geneid"].isin(Geneids_count0))
            &
            (new_coordinates["score_peakmean"]>0)
            &
            (new_coordinates["dist_to_TandemGene"]<200) #next gene is less than 200bp away
            &
            (
                ((new_coordinates["peak_to_start"]<-100) & (new_coordinates["strand"]=="-"))
                |
                ((new_coordinates["peak_to_end"]>100) & (new_coordinates["strand"]=="+"))
                )
                )
        
        # sets the value of novel transcript to True if CDS_count0 is True
        # new_coordinates["novel_transcript"] = new_coordinates.apply(lambda row: \
        #     True if row["CDS_count0"]\
        #         else row["novel_transcript"],
        #         axis = 1
        #     )

        return new_coordinates


    def hypothetical_correction(self, new_coordinates):
        """
        adds row for genes without peak, that were initially assigned hypothetical
        does not overwrite hypothetical annotation

        """

        # some genes were overwritten as hypothetical, this code adds a new line for each overwritten gene
        hypothetical_only = new_coordinates[(new_coordinates["new_attribute"].str.contains("h")) & (~new_coordinates["attributevalue"].duplicated())]
        hypothetical_only[
            ["novel_transcript","no_peak", "score_peak",
             "score_peakmean", "peak_to_start", "peak_to_end",
                # "blacklist"
                ]
            ] = False,True, -20, -20, 0, 0
        
        hypothetical_only["peak_locationmin"] = hypothetical_only["start"]
        hypothetical_only["peak_locationmax"]= hypothetical_only["end"]
        hypothetical_only["new_attribute"] = hypothetical_only["new_attribute"].str.replace("h", "")#.str.replace("bl", "")
        new_coordinates = new_coordinates.append(hypothetical_only)
        
        
        return new_coordinates

    
    def plot_result(self):
        
        if not os.path.isfile(self.PS_data_output):
            new_coordinates = self.CollapseGTF()
        else:
            try:
                new_coordinates = self.new_coordinates
            except Exception as e:
                new_coordinates = pd.read_csv(self.PS_data_output)

        fig, axes =plt.subplots(figsize = (10,5), ncols = 2, sharey= True)
        markers = {True: "X", False: "o"}
        dict_strand = dict(zip(["-", "+"], [-1, 1]))
        
        for ax, peak_to, strand in zip(axes, ["start", "end"], ["-", "+"]):
            data = new_coordinates[new_coordinates["strand"] ==strand].reset_index()
            # data["blacklist_actual"] = (data["blacklist"]) & (abs(data[f"peak_to_{peak_to}"]) > 500)


            sns.scatterplot(data = data, x =f"peak_to_{peak_to}", y = "score_peakmean", ax =ax, style = "novel_transcript",
                            # hue = "blacklist_actual", 
                            # hue_order = [False, True],  
                            # palette = ["#012536", "#f04175"],
                            color = "#012536",
                            markers =markers,
                            s= 20, legend = True, alpha = .5)
            ax.set_yscale("log")
            for x in [1000,500,-400]:
                x= dict_strand[strand]*x
                ax.vlines(ymin = 1, ymax = 10**7, x = x, color = "k", lw = 0.5)
                ax.text(x = x, y = 10**7, s = f"{x}", fontsize=6)
            ax.hlines(xmin = ax.get_xlim()[0], xmax = ax.get_xlim()[1], y = np.mean(self.df_summary["q75"]), color = "k", lw = 0.5)
            dict_strand = dict(zip(["-", "+"], [-1, +1]))

            ax.set_xlim(np.sort((dict_strand[strand]*5000, dict_strand[strand]*-1000)))
            ax.set_ylabel("peak score (mean)")
            ax.set_xlabel(f"distance of ({strand}) peak to {peak_to} of ({strand}) gene")

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

        if os.path.isfile(self.new_coordinates_file_path) and not self.rerun:
            raise FileExistsError(f"Output file already exists: {self.new_coordinates_file_path}")
        
        print("Finalizing assignments of peaks to gene coordinates...")
        self.RunAll()

        print("Matching peaks to gene annotations...")
        peak_links2GTF = self.MatchtoGTF()
        peak_links2GTF = self.Attribute2GTF(peak_links2GTF)
        new_coordinates = self.aggregatorGTF(peak_links2GTF)

        print("Corrrecting counts based on CDS counts, and reasssigning attribute...")
        new_coordinates = self.CDScounts_correction(new_coordinates)
        new_coordinates = self.Attribute2GTF(new_coordinates, attribute_col = "attributevalue")
        new_coordinates = self.hypothetical_correction(new_coordinates)

        new_coordinates.sort_values(["seqname", "start", "end"], inplace=True)
        new_coordinates["gene_size"] = abs(new_coordinates["start"] -  new_coordinates["end"])
        
        # use specified gene_biotypes (e.g. protein_coding for new annotation file)
        new_coordinates = new_coordinates[new_coordinates["gene_biotype"].isin(self.gene_biotypes)]
        
        new_coordinates.rename(columns = {"new_attribute":"attribute"}, inplace=True)
        new_coordinates.drop_duplicates(["seqname", "strand", "attribute"], keep = "last", inplace = True)
        new_coordinates["new_start"]  = new_coordinates.apply(lambda row: self.new_start(row) , axis =1)
        new_coordinates["new_end"]  = new_coordinates.apply(lambda row: self.new_end(row) , axis =1)
        new_coordinates["start_diff"] = new_coordinates["new_start"] -  new_coordinates["start"]
        new_coordinates["end_diff"] = new_coordinates["new_end"] -  new_coordinates["end"]

        

        new_coordinates = self.MakeSaveGTF(new_coordinates)
        self.SaveNewGTF(new_coordinates)
        
        return new_coordinates 


    def MakeSaveGTF(self, new_coordinates):
        print("Making new coordinates file...")


        
        new_coordinates["start"] = new_coordinates["new_start"].astype(int)
        new_coordinates["end"] = new_coordinates["new_end"].astype(int)
        new_coordinates["frame"] = "."
        new_coordinates["feature"] = "predicted_3primeUTR"
        new_coordinates["source"] = "PeakStream"

        print("Annotating new coordinates...")
        new_coordinates.to_csv(self.PS_data_output, index = False)
        print(f"Saving new coordinates data to {self.PS_data_output}")
        

    
        return new_coordinates
    
   
    def RunAll(self):
        self.putativeUTRS()
        self.AT_stretch()
        self.countpeaks()
        self.aggregateLinks() 
        self.applySummary()
        self.setRules()
        self.applyRules()
        self.aggregateThresh()
        self.Collapse()

        return self.peak_links

    
    def SaveNewGTF(self, new_coordinates):

        GTF_cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand','frame', 'attribute']
        new_coordinates_save = new_coordinates[GTF_cols]

        new_coordinates_save.to_csv(self.new_coordinates_file_path, sep = "\t", index=False, header = None,  quoting=csv.QUOTE_NONE)
        print(f"PeakStream run completed. \nSaved new GTF file to {self.new_coordinates_file_path}")
        
        return self.new_coordinates_file_path


    def PeakStream2featureCounts(self):


        if not os.path.isfile(self.new_coordinates_file_path) or self.rerun:
            self.CollapseGTF()
            
        
        User_inputs_dict_update = self.User_inputs_dict
        User_inputs_dict_update["genome_annotations_path"] = self.new_coordinates_file_path 
        featurecounts_file = featurecounts(self.directory, User_inputs_dict_update, s=1, t="predicted_3primeUTR", T=20, 
                                        stardir="starout", featurecounts_subDir="featurecounts", rerun=self.rerun)
        featurecounts_file_biotype = featurecounts_biotype(featurecounts_file, self.gene_biotypes, t="gene", PeakStream =True)
        

        return featurecounts_file_biotype
            