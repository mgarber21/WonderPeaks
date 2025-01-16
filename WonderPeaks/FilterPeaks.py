"""
The following functions are applied to the WonderPeaks 
unfiltered output to filter for peaks that meet 
relevant conditions

"""


import shutil
import sys
import os


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from collections import defaultdict



# from scipy.stats import hypergeom
# from scipy.misc import electrocardiogram
# from scipy.signal import find_peaks
# from scipy.signal import argrelextrema

# import math
from itertools import combinations, permutations, groupby

from user_inputs import *
from Helper_Functions import *
from WonderPeaks.WonderPeaks4ChIP import *


class Filter_WonderPeaks(WonderPeaks4ChIP):
    """
    Final processing of PeakStream    
    
    """
    
    def __init__(self, directory, coordinate_file, 
                       designfactors, slope_thresh,
                 ):
        
        
        
        
        self.directory = directory
        self.WP_directory = WonderPeaks_Direcory(directory)
        self.coordinate_file = coordinate_file
        self.designfactors = designfactors
        self.slope_thresh =slope_thresh
        
        self.WOnder_init = os.path.join(self.WP_directory, "WOnder_init.csv")
        self.WOnder_unfiltered_peaks = os.path.join(self.WP_directory, "WOnder_unfiltered_peaks.csv")
        self.WPfigureDir = os.path.join(self.WP_directory, "WonderPeaks_figures")
        os.makedirs(self.WPfigureDir, exist_ok=True)


        WonderPeaks4ChIP.__init__(self, directory, coordinate_file, designfactors, slope_thresh)
        self.WPC = WonderPeaks4ChIP(directory, coordinate_file,designfactors,slope_thresh)
        self.df_peaks= self.call_peaks()
        # self.init_data = pd.read_csv(self.WOnder_init)




    def metadata_connect(self):
        metadata = self.metadata
        metadata["replicate"] = metadata.groupby("_".join(self.designfactors)).cumcount() + 1
        self.df_peaks  = pd.merge(self.df_peaks, self.metadata, on = "basename")
        
        return self.df_peaks

    def get_tag_factor(self):
        designfactors_col = "_".join(self.designfactors)

        designfactor_dict =  dict()
        for designfactor in self.designfactors:
            factor_samples = self.metadata[designfactor].unique()
            if [factor_sample for factor_sample in factor_samples if "tag" in factor_sample.lower()]:
                # tagged factor is designfactor column to use for lowest level
                tagged_factor = designfactor
            else:
                other_factor = designfactor
            designfactor_dict[designfactor] = list(factor_samples)
        return designfactor_dict, tagged_factor, other_factor

    def make_tagged_factor_rep(self):
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()
        self.df_peaks[f"{tagged_factor}_rep"] = self.df_peaks.apply(lambda row: row[tagged_factor]+"_"+str(row["replicate"]),
                                        axis= 1)
        return self.df_peaks
    
    


    def get_sample_combos(self, data):
        
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()
        ## create dictionaries using sample metadata
        #  Note control will always be first in the dictionart because dictionaries are sorted alphabetically (control, test)
        dict_sample_types = {}
        dict_sample_types["control"] = [
            i for i in data[tagged_factor].unique() if "control" in i.lower() or "ctrl" in i.lower()
            ][0]
        dict_sample_types["test"]= [
            i for i in data[tagged_factor].unique() if not ("control" in i.lower() or "ctrl" in i.lower())
            ][0]
        
        dict_sample_types_reps = {}
        dict_sample_types_reps["control"] = [i for i in data[f"{tagged_factor}_rep"].unique()if "control" in i.lower() or "ctrl" in i.lower() or "untagged" in i.lower() or "input" in i.lower()]
        dict_sample_types_reps["test"] = [i for i in data[f"{tagged_factor}_rep"].unique() if not ("control" in i.lower() or "ctrl" in i.lower() or "untagged" in i.lower() or "input" in i.lower())]
        
        
        #  Note control will always be first because dictionaries are sorted
        

      
        #  Note control will always be first because dictionaries are sorted
        Sample_type_combos = list(combinations(data[f"{tagged_factor}_rep"].unique(), 2))
        
        dict_sample_type_combos = {}
        for sample_type in list(dict_sample_types.values()) + [list(dict_sample_types.values())] :
            if type(sample_type) == list:
                combos_return = [i if (i[0].startswith(sample_type[0]) and i[1].startswith(sample_type[1]))
                                 else i[::-1]
                                 for i in Sample_type_combos
                                 if (i[0].startswith(sample_type[0]) and i[1].startswith(sample_type[1]))
                                 or (i[0].startswith(sample_type[1]) and i[1].startswith(sample_type[0]))
                                 ]
                dict_sample_type_combos["V".join(sample_type)] = combos_return
                
                # dict_sample_type_combos["V".join(sample_type)] = [i for i in Sample_type_combos if i[0].startswith(sample_type[0]) and i[1].startswith(sample_type[1])]
            else:    
                dict_sample_type_combos[f"{sample_type}V{sample_type}"] = [i for i in Sample_type_combos if i[0].startswith(sample_type) and i[1].startswith(sample_type)]

        return dict_sample_types,dict_sample_types_reps, dict_sample_type_combos


    def other_factor4file_name(self, other_factor):
        other_factor_value = data[other_factor].unique()
        return other_factor_value

    def plot_raw_results(self, save = True, figsize_x = 10, figsize_y=3):
        
        fig, ax = plt.subplots(figsize = (figsize_x,figsize_y))
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()

        ax = sns.stripplot(data = self.df_peaks, y = "normalized_score_peak", x = f"{tagged_factor}_rep", 
                    hue = other_factor,
                    alpha = 0.25, dodge=True, s = 2.5,
                    palette= "Greys",
                    ax=ax)



        ax.set_yscale("log")


        ax.grid(visible =True, axis = "y", alpha = 0.3)

        ax.legend(bbox_to_anchor = (1,1))
        plt.draw()

        xticklabels = ax.get_xticklabels()
        new_labels = ["\nreplicate".join(label.get_text().replace("_", " ").rsplit(" ", 1)) for label in xticklabels]
        ax.set_xticks(range(len(new_labels)))
        ax.set_xticklabels(new_labels)
        ax.set_xlabel("")
        ax.set_ylabel("Peak Score\n(normalized to average score)")
        ax.set_title("Raw Data Summary")
        plt.tight_layout()
        if save:
            
            plt.savefig(os.path.join(self.WPfigureDir, "raw_data_summary.png"))
        else:
            plt.show()


    def filter_peaks(self, 
                    designfactor_value, 
                    score_cut,
                    fold_greater, 
                    ):
        self.designfactor_value = designfactor_value
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()
        dict_sample_types,  dict_sample_types_reps, dict_sample_type_combos = self.get_sample_combos(self.df_peaks)

        data_value = self.df_peaks[self.df_peaks[other_factor] == self.designfactor_value].reset_index().copy()
        data_value.sort_values(["seqname", "peak_location", "normalized_score_peak"], inplace=True)


        dict_sample_types,  dict_sample_types_reps, dict_sample_type_combos = self.get_sample_combos(data_value)


        data_value["duplicate_index"] = data_value.groupby([other_factor, "seqname", "peak_location", f"{tagged_factor}_rep",])["peak_location"].cumcount()

        data_filter = data_value.pivot(
            columns = tagged_factor + "_rep", 
            index= [other_factor, "seqname", "peak_location", "duplicate_index"],
            values = "normalized_score_peak").fillna(1)
        
        data_filter.reset_index(inplace=True)

        data_filter["median_control"] = data_filter.apply(
                    lambda row: np.median(row[dict_sample_types_reps["control"]]),
                    axis=1)
        data_filter["median_test"] = data_filter.apply(
                    lambda row: np.median(row[dict_sample_types_reps["test"]]),
                    axis=1)

        
        data_filter["fold_change"] = data_filter["median_test"]/data_filter["median_control"]
        data_filter["log2fold_change"] = np.log2(data_filter["fold_change"])
                    
        
        # # Add a new column based on the condition
        data_filter["unique_peak"] = data_filter.apply(
            lambda row: (
                
                
                True if (
                    all(row[test]/row[control] > fold_greater\
                        for control, test in dict_sample_type_combos[f"{dict_sample_types.get('control')}V{dict_sample_types.get('test')}"]
                        ) and
                    any(row[control] > 1 and row[test] > score_cut \
                        for control, test in  dict_sample_type_combos[f"{dict_sample_types.get('control')}V{dict_sample_types.get('test')}"]
                        )
                )
                else True if (
                    ((row["fold_change"] > fold_greater) and (len(dict_sample_types_reps["test"])> 2)) and 
                    all(row[control] >= 1 and row[test] > score_cut \
                        for control, test in  dict_sample_type_combos[f"{dict_sample_types.get('control')}V{dict_sample_types.get('test')}"]
                        )
                )
                else True if (
                    all((row[control]== 1 and row[test] > 3)\
                        for control, test in  dict_sample_type_combos[f"{dict_sample_types.get('control')}V{dict_sample_types.get('test')}"]
                        )
                )

                else False if (
                    any(row[control] > 1 and row[test] > score_cut\
                        for control, test in  dict_sample_type_combos[f"{dict_sample_types.get('control')}V{dict_sample_types.get('test')}"]
                        )

                ) 
                
                else True if (
                    any(row[test_i] > score_cut and row[test_j] > score_cut \
                        for test_i, test_j in dict_sample_type_combos[f"{dict_sample_types.get('test')}V{dict_sample_types.get('test')}"]) 
            )
                            
                else False
            ),
            axis=1
        )
        
        
            
        return data_value, data_filter

    def plot_tVut(self, data_filter, save = True):
        
            dict_sample_types,  dict_sample_types_reps, dict_sample_type_combos = self.get_sample_combos(self.df_peaks)
            # Grouping by the first element of each tuple
            grouped_data = defaultdict(list)
            for x in sorted(dict_sample_type_combos[f"{dict_sample_types.get('control')}V{dict_sample_types.get('test')}"]):
                grouped_data[x[0]].append(x)

            # Convert to a list of lists
            grouped_data = list(grouped_data.values())
            # grouped_data = [list(group) for key, group in groupby(dict_sample_type_combos[f"{dict_sample_types.get('control')}V{dict_sample_types.get('test')}"] ,
                                                            # key=lambda x: x[0])]
            # return grouped_data
            fig, axes = plt.subplots(figsize = (3*len(grouped_data[0]),3*len(grouped_data)) , ncols = len(grouped_data), nrows = len(grouped_data[0]))

            for xys, ax_rows in zip(grouped_data, axes):
                for xy, ax in zip(xys, ax_rows): 
                        
                        data_plot = data_filter[data_filter[xy[0]]+data_filter[xy[1]]!=0]
                        data_plot = data_plot[~((data_plot[xy[0]]==1) & (data_plot[xy[1]]==1))]
                        
                        sns.scatterplot(
                        data = data_plot, 
                        y =xy[0], x = xy[1],
                        style = "unique_peak", style_order = [True, False],
                        hue = "unique_peak", hue_order = [True, False],
                        palette = ["#f0cf81", "black"], alpha = .5 ,s=12,
                        ax =ax, 
                        legend = False, 
                                        )
                        
                        data_plot2 =data_plot[data_plot["unique_peak"]]
                        # bring unique peaks to front of figure
                        sns.scatterplot(
                        data = data_plot2, 
                        y =xy[0], x = xy[1],
                        color = "#f0cf81", alpha = .15 ,s=12,
                        ax =ax, edgecolor = "#f59d40",
                        legend = False, 
                                        )
                        ax.set_xscale("log")
                        ax.set_yscale("log")
                        lim = max(ax.get_ylim()[1], ax.get_xlim()[1])
                        ax.set_ylim(0.9, lim)
                        ax.set_xlim(0.9, lim)
                        ax.plot([0, 1], [0, 1], transform=ax.transAxes, color = "black", lw = .5, 
                                alpha = 0.5, zorder = 0)
                        xlabel, ylabel = (
                                "\n replicate ".join(ax.get_xlabel().rsplit("_", 1)),
                                "\n replicate ".join(ax.get_ylabel().rsplit("_", 1)),
                                )
                        xlabel, ylabel = [" ".join(label.split("_")) for label in (xlabel, ylabel)]
                                
                        ax.set_ylabel(ylabel, fontsize=10)
                        ax.set_xlabel(xlabel, fontsize=10)

            plt.tight_layout()
            if save:
                plt.savefig(os.path.join(self.WPfigureDir, f"{self.designfactor_value}_untaggedVtagged.png"))
                # plt.savefig(os.path.join(directory, "publication_figures", f"Supplemental_allVall_{other_factor_value}_untaggedVtagged.pdf"), format ="pdf")
            plt.show()


    def plot_tVt_and_utVut(self, data_filter, save = True):
            dict_sample_types,  dict_sample_types_reps, dict_sample_type_combos = self.get_sample_combos(self.df_peaks)
            control_grouped_data = dict_sample_type_combos[f"{dict_sample_types.get('control')}V{dict_sample_types.get('control')}"] 
            test_grouped_data = dict_sample_type_combos[f"{dict_sample_types.get('test')}V{dict_sample_types.get('test')}"] 

            fig, axes = plt.subplots(figsize = (3*len(test_grouped_data),6) , ncols = len(test_grouped_data), nrows = 2)

            for xys, row_ax in zip([control_grouped_data, test_grouped_data], axes):
                    for xy, ax in zip(xys, row_ax):
                            data_plot = data_filter[data_filter[xy[0]]+data_filter[xy[1]]!=0]
                            data_plot = data_plot[~((data_plot[xy[0]]==1) & (data_plot[xy[1]]==1))]
                            
                            sns.scatterplot(
                            data = data_plot, 
                            y =xy[0], x = xy[1],
                            style = "unique_peak", style_order = [True, False],
                            hue = "unique_peak", hue_order = [True, False],
                            palette = ["#fab12a", "black"], alpha = .5 ,s=12,
                            
                            ax =ax, 
                            legend = False,
                            )
                            data_plot2 =data_plot[data_plot["unique_peak"]]
                            # bring unique peaks to front of figure
                            sns.scatterplot(
                            data = data_plot2, 
                            y =xy[0], x = xy[1],
                            color = "#f0cf81", alpha = .15 ,s=12,
                            ax =ax, edgecolor = "#f59d40",
                            legend = False, 
                                            )
                            
                            ax.set_xscale("log")
                            ax.set_yscale("log")
                            lim = max(ax.get_ylim()[1], ax.get_xlim()[1])
                            ax.set_ylim(0.9, lim)
                            ax.set_xlim(0.9, lim)
                            ax.plot([0, 1], [0, 1], transform=ax.transAxes, color = "black", lw = .5, 
                            alpha = 0.5, zorder = 0)
                            xlabel, ylabel = (
                                    "\n replicate ".join(ax.get_xlabel().rsplit("_", 1)),
                                    "\n replicate ".join(ax.get_ylabel().rsplit("_", 1)),
                                    )
                            xlabel, ylabel = [" ".join(label.split("_")) for label in (xlabel, ylabel)]
                                    
                            ax.set_ylabel(ylabel, fontsize=10)
                            ax.set_xlabel(xlabel, fontsize=10)
                            
                    
            plt.tight_layout()
            if save:
                plt.savefig(os.path.join(self.WPfigureDir, f"{self.designfactor_value}_taggedVtagged_untaggedVuntagged.png"))
                # plt.savefig(os.path.join(directory, "publication_figures", f"Supplemental_allVall_{other_factor_value}_untaggedVtagged_taggedVtagged.pdf"), format ="pdf")

            plt.show()
    def get_baseline(self, file, seqname, peak_location, peak_diff, median):
        """Extracts baseline score from bedgraph file of nearby point"""
        df_bed = pd.read_csv(file, sep = "\t")
        df_bed.rename(columns=dict(zip(df_bed.columns, ["seqname", "start", "stop", "score"])), inplace =True)
        bed_convert_up,bed_convert_down  = np.round(((peak_location-peak_diff)/2)/20)*20, np.round(((peak_location+peak_diff+20)/2)/20)*20
        df_bed_loc_up = df_bed[(df_bed["seqname"] == seqname) &(df_bed["start"]== max(20, bed_convert_up))]
        df_bed_loc_down = df_bed[(df_bed["seqname"] == seqname) &(df_bed["start"]== max(20, bed_convert_down))]

        if df_bed_loc_up.empty and df_bed_loc_down.empty:
            score_up, score_down = median, median
            return score_up, score_down
        elif df_bed_loc_up.empty:
            score_up, score_down = median, df_bed_loc_down["score"].iloc[0]
            return score_up, score_down
        elif df_bed_loc_down.empty:
            score_up, score_down = df_bed_loc_up["score"].iloc[0], median
            return score_up, score_down
        else:
            score_up, score_down = df_bed_loc_up["score"].iloc[0], df_bed_loc_down["score"].iloc[0]
            return score_up, score_down
        
        
    def write_baseline(self, peaks):
        """Applys get baseline to rows of peaks dataframe"""
        baselines = peaks.apply(lambda row: 
            self.get_baseline(
               row["file"], row["seqname"], row["peak_location"], row["peak_diff"], row["median"])
            ,axis=1
        )
        return baselines
    
    def label_false_positives(self, peaks, tagged_factor, enrichment_thresh):
        """
        Applies functions to label or remove peaks from noisy genomic regions. 
        False positives are multiple consecutive high peaks

        returns: 
        - pd.DataFrame: DataFrame of peak coordinates in your bedgraph files above a certain width threshold
        
        """

        peaks["peak_location_group_diff"] = peaks.groupby(["seqname", "peak_location"])["peak_location_actual"].transform(lambda x: x.max()-x.min()).astype(int)
        peaks.fillna(1, inplace=True)
        
        # peaks["baselines"] = self.write_baseline(peaks)
        peaks["label_false_pos"] = peaks.apply(lambda row: 
           True if (row["peak_location_group_diff"]>2*row["peak_diff"])
           else True if (abs(row["1st_derivative"])<0.2 and  row["score_peak"]/max(1, row["score_left_edge"]) < enrichment_thresh and row["score_peak"]/max(1, row["score_left_edge"]) < enrichment_thresh)
           else True if (row["peak_diff"]>300 and abs(row["1st_derivative"])<0.3 and  row["score_peak"]/max(1, row["score_left_edge"]) < enrichment_thresh and row["score_peak"]/max(1, row["score_left_edge"]) < enrichment_thresh)
        #    else True if (row["score_peak"]-row["baselines"][0] < 10  and row["score_peak"]-row["baselines"][1] < 10)
           else False,
           axis = 1)
        
        # anything that got labeled as True should extend to the group
        peaks['label_false_pos'] = peaks.groupby(["seqname", "peak_location", tagged_factor])['label_false_pos'].transform(lambda x: x.any())

        return peaks
    
    def data_TruePeaks(self, data_value, data_filter,enrichment_thresh ):
        
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()
        dict_sample_types,  dict_sample_types_reps, dict_sample_type_combos = self.get_sample_combos(data_value)
        

        TruePeaks_piv = data_filter[data_filter["unique_peak"]]
        
        TruePeaks = TruePeaks_piv.melt(
                id_vars=['Cell_type', 'seqname', 'peak_location', "fold_change"],
                value_vars=dict_sample_types_reps["test"],
                value_name = "normalized_score_peak")
        
        TruePeaks["log2fold_change"] = np.log2(TruePeaks["fold_change"])
        TruePeaks[tagged_factor] = TruePeaks.apply(
                lambda row: row[f"{tagged_factor}_rep"].split("_")[0],
                axis=1)

        for method in ["mean", "median", "min", "max", "std" ]:
                TruePeaks[method] =  TruePeaks.groupby(['seqname', 'peak_location', tagged_factor])["normalized_score_peak"].transform(method)


        TruePeaks = TruePeaks[(TruePeaks["median"] > 2) |
                                    (TruePeaks["min"] != 1)
                                    ].reset_index()

        TruePeaks = pd.merge(
                        data_value[["file","seqname", "median","peak_location", "peak_location_actual", 
                                    other_factor,tagged_factor, f"{tagged_factor}_rep","score_peak",
                                    "normalized_score_peak", "peak_location_broad",
                                    "1st_derivative", "score_left_edge", "score_right_edge",
                                    "peak_diff", "start_peak"]],
                        TruePeaks,
                        on = ["seqname", other_factor, tagged_factor, f"{tagged_factor}_rep","peak_location", "normalized_score_peak"],
                        suffixes = ("", "True")
                        )
        TruePeaks = self.label_false_positives(TruePeaks, tagged_factor, enrichment_thresh)
        TruePeaks["peak_location_start"] = TruePeaks["peak_location_actual"] - TruePeaks["peak_diff"] 
        TruePeaks["peak_location_stop"] = TruePeaks["peak_location_actual"] + TruePeaks["peak_diff"] +20 
        
        return TruePeaks




    def load_coordinates_withnames(self):
        gtf =load_coordinates(self.coordinate_file, intergenic = True, max_intergenic_dist= 2000)
        for named_attribute in ["gene_id", "gene_name"]:
            gtf[named_attribute] = gtf.apply(
                lambda row: rename_attributes(row["attribute"], named_attribute = named_attribute)\
                    if rename_attributes(row["attribute"], named_attribute = named_attribute)
                    else rename_attributes(row["attribute"], named_attribute = "gene_id")
                
                ,
                axis = 1
            )
        return gtf

    def Peaks2Annotations(self, data_True):
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()

        gtf=self.load_coordinates_withnames()
        peaks2gtf = pd.DataFrame()

        for  method, right_on  in zip(["forward", "backward"], ["intergenic_end", "intergenic_start"]):
            # gtf_strand = gtf[gtf["strand"] == strand]
            gtf[right_on] = gtf[right_on].astype("int")
            gtf.sort_values(by =[right_on], inplace=True)
            data_True["peak_location"] = data_True["peak_location"].astype("int")
            data_True.sort_values(by =["peak_location"], inplace=True)


            data_True_temp = pd.merge_asof(data_True, gtf, 
                        left_index=False, right_index=False, 
                        left_on="peak_location", right_on=right_on, 
                        by=["seqname"],
                        suffixes = ("_peaks", ""), tolerance=20000, 
                        allow_exact_matches=True, direction=method)
            
            peaks2gtf = pd.concat([peaks2gtf, data_True_temp])
        peaks2gtf["peak_in_intergenic"] = peaks2gtf.apply(lambda row:
                                            True if row['peak_location'] >= row["intergenic_start"] and row['peak_location'] <= row["intergenic_end"]\
                                            else False,
                                            axis =1
                                            )
        peaks2gtf = peaks2gtf[peaks2gtf["peak_in_intergenic"]].reset_index(drop=True)
        peaks2gtf.sort_values(["seqname", "peak_location"], inplace=True)
        peaks2gtf["peaks_per_gene"] = peaks2gtf.groupby(["gene_name",f"{tagged_factor}_rep" ])["peak_location"].transform("count")
        peaks2gtf["peak_distance_to_gene"] = peaks2gtf.apply(lambda row:  
            row["intergenic_end"] - row["peak_location"] \
                if (row["strand"] == "+") \
                    else -(row["peak_location"]-row["intergenic_start"]),
                    axis =1
                    )
        return peaks2gtf


    def summary_dist_to_gene(self, peaks2gtf, save = True):
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()
        ncols = len(peaks2gtf[f"{tagged_factor}_rep"].unique())
        fig, axes = plt.subplots(figsize = (4*ncols,4), ncols = ncols)
        for ax, rep in zip(axes, sorted(peaks2gtf[f"{tagged_factor}_rep"].unique())):
            plot_data = peaks2gtf[
                    (peaks2gtf[f"{tagged_factor}_rep"] == rep) & (~peaks2gtf["label_false_pos"])]
            
            sns.histplot(data =plot_data , 
                        x = "peak_distance_to_gene",
                        hue = "strand", 
                        palette= ["#a298fa", "#fa98ce"],
                        ax =ax, kde =True,element="step"
                        )
            ax.set_title(" ".join(rep.split("_")))
            ax.set_xlabel(" ".join(ax.get_xlabel().split("_")))
        fig.suptitle(f"{self.designfactor_value}\nDistance of peaks to TSS")
        # fig.subplots_adjust(top=0.85)
        plt.tight_layout()
        if save:
            plt.savefig(os.path.join(self.WPfigureDir, f"{self.designfactor_value}_peak_dist2gene_summary.png"))
            
        

    def summary_peaks_per_gene(self, peaks2gtf, text_thresh = 20, save= True, log2 = True):
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()

        data = peaks2gtf.drop_duplicates(subset = ["peak_location", "gene_name"], keep="first")
        data.reset_index(drop =True, inplace=True)
        

        x = "peaks_per_gene"
        y= "fold_change"
        if log2:
            y="log2fold_change"
            text_thresh = np.log2(text_thresh)
        
        
        x_max = np.max(data[x])


        fig, ax = plt.subplots(figsize = (max(4, x_max),5) )

        sns.scatterplot(data = data, x = x, y=y ,
                        color = "white",
                        edgecolor="black",
                
                        alpha = .5 ,s=25,
                        ax =ax, 
                        legend = False,)

        ax.hlines(xmin =ax.get_xlim()[0], xmax=ax.get_xlim()[1], y = text_thresh, zorder = 0, linestyles="--", color = "grey", lw = 1)

        
        text = data[(data[y]> text_thresh) & (~data["label_false_pos"])]\
                .reset_index(drop =True)\
                        .sort_values([y])

        sns.scatterplot(data = text, x = x, y= y,
                        color = "#5363db",
                        edgecolor="black",
                        alpha = .5 ,s=25,
                        ax =ax, 
                        legend = False, zorder=1)

        text_range = int(ax.get_xlim()[1]/50)
        x0, y0 = 0,0
        texts = []
        for _, row in text[(text["strand"] == "+")].iterrows():
                x1, y1 = row[x], row[y]
                if ((abs(y0-y1) < ax.get_ylim()[1]*0.5) and (abs(x0-x1) < 1)):
                        x0, y0 = x1, y1
                        continue
                x0, y0 = x1, y1
                x2, y2 = row[x]+text_range, row[y]
                annotation = row["gene_name"]
                texts.append(ax.annotate(annotation,# Annotation text
                        xy = (x1, y1),  # x,y-coordinate
                        xytext=(x2, y2), #x,y-offset
                        fontsize=8,  # Font size
                        ha='left',  # Horizontal alignment
                        va='bottom',  # Vertical alignment
                        color='#fa98ce',  # Text color (optional),
                        alpha=1,
                                ))
        x0, y0 = 0,0
        for _, row in text[(text["strand"] == "-")].iterrows():
                x1, y1 = row[x], row[y]
                if ((abs(y0-y1) < 4) and (abs(x0-x1) < 1)):
                    x0, y0 = x1, y1
                    continue
                x0, y0 = x1, y1
                x2, y2 = row[x], row[y]
                annotation = row["gene_name"]
                texts.append(ax.annotate(annotation,# Annotation text
                        xy = (x1, y1),  # x,y-coordinate
                        xytext=(x2, y2), #x,y-offset
                        fontsize=8,  # Font size
                        ha='right',  # Horizontal alignment
                        va='top',  # Vertical alignment
                        color='#a298fa',  # Text color (optional),
                        alpha=1,
                                ))
        
        false_pos = data[data["label_false_pos"]].reset_index(drop =True)

        sns.scatterplot(data = false_pos, x = x, y= y,
                        color = "#cf063c",
                        edgecolor="black",
                        alpha = .5 ,s=25,
                        ax =ax, 
                        legend = False,)
        
        # Get unique x values and sort them
        unique_x = sorted(data[x].unique())

        # Set the x-axis ticks to unique values
        ax.set_xticks(unique_x)

        ax.set_title("Number of peaks associated with a gene\nvs.\nPeak Fold Change")
        ax.set_xlabel(" ".join(ax.get_xlabel().split("_")))
        ax.set_ylabel(" ".join(ax.get_ylabel().split("_")))
        ax.spines[['right', 'top']].set_visible(False)

        ax.margins(x=0.05)
        fig.suptitle(f"{self.designfactor_value}", ha='left')    
        plt.tight_layout()
    
        if save:
                plt.savefig(os.path.join(self.WPfigureDir, f"{self.designfactor_value}_peaks_per_gene_summary.png"))
        plt.show()



    def summary_ranks(self, peaks2gtf, text_thresh = 25, save= True, log2 = True):
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()

        data = peaks2gtf.copy().drop_duplicates(["peak_location", "gene_name"])
        data.reset_index(drop =True, inplace=True)
        data["rank"] = data["fold_change"].rank(ascending=False, method="min")

        fig, ax = plt.subplots(figsize = (10,5) )
        x = "rank"
        y= "fold_change"
        if log2:
            y="log2fold_change"
            text_thresh = np.log2(text_thresh)

        sns.scatterplot(data = data, x = x, y=y ,
                        color = "white",
                        edgecolor="black",
                
                        alpha = .5 ,s=25,
                        ax =ax, 
                        legend = False,)

        ax.hlines(xmin =ax.get_xlim()[0]-0.5, xmax=ax.get_xlim()[1], y = text_thresh, zorder = 0, linestyles="--", color = "grey", lw = 1)



        text = data[(data[y]> text_thresh) & (~data["label_false_pos"])].reset_index(drop =True)
        text.sort_values([x, y], inplace=True)

        sns.scatterplot(data = text, x = x, y= y,
                        color = "#5363db",
                        edgecolor="black",
                        alpha = .5 ,s=25,
                        ax =ax, 
                        legend = False,
                        zorder=1)

        x0, y0 = 0,0
        text_range = int(ax.get_xlim()[1]/50)
        texts = []

        for _, row in text[(text["strand"] == "+")]\
                .reset_index(drop=True)\
                .iterrows():
                x1, y1 = row[x], row[y]
                if (x0, y0 == (0,0)):
                    pass
                if (abs(y0-y1) < ax.get_ylim()[1]*0.5 and abs(x0-x1) < ax.get_ylim()[1]*0.05):
                        x0, y0 = x1, y1
                        continue
                x0, y0 = x1, y1
                x2, y2 = row[x]+text_range, row[y]
                annotation = row["gene_name"]
                texts.append(ax.annotate(annotation,# Annotation text
                        xy = (x1, y1),  # y-coordinate
                        xytext=(x2+0.5, y2),
                        fontsize=8,  # Font size
                        ha='left',  # Horizontal alignment
                        va='bottom',  # Vertical alignment
                        color='#fa98ce',  # Text color (optional),
                        alpha=1,
                                ))
        x0, y0 = 0,0
        for _, row in text[(text["strand"] == "-")]\
                .reset_index(drop=True)\
                        .iterrows():
                x1, y1 = row[x], row[y]
                if abs(y0-y1) < 4 and abs(x0-x1) < 3:
                        x0, y0 = x1, y1
                        continue
                x0, y0 = x1, y1
                x2, y2 = row[x]-text_range, row[y]
                annotation = row["gene_name"]
                texts.append(ax.annotate(annotation,# Annotation text
                        xy = (x1, y1),  # y-coordinate
                        xytext=(x2-0.5, y2),
                        fontsize=8,  # Font size
                        ha='right',  # Horizontal alignment
                        va='bottom',  # Vertical alignment
                        color='#a298fa',  # Text color (optional),
                        alpha=1,
                                ))

                
        false_pos = data[data["label_false_pos"]].reset_index(drop =True)

        sns.scatterplot(data = false_pos, x = x, y= y,
                        color = "#cf063c",
                        edgecolor="black",
                        alpha = .5 ,s=25,
                        ax =ax, 
                        legend = False,)
        ax.set_title("Gene Rank vs. Peak Fold Change")
        ax.set_xlabel(" ".join(ax.get_xlabel().split("_")))
        ax.set_ylabel(" ".join(ax.get_ylabel().split("_")))
        ax.spines[['right', 'top']].set_visible(False)

        ax.margins(x=0.05)
        fig.suptitle(f"{self.designfactor_value}", ha='left')
        plt.tight_layout()       
        if save:
                plt.savefig(os.path.join(self.WPfigureDir, f"{self.designfactor_value}_peak_ranks_summary.png"))
        plt.show() 
            
            
    def peak_locations2narrowPeak(self, data):
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()
        
        for rep in data[f"{tagged_factor}_rep"].unique():
            data_rep = data[data[f"{tagged_factor}_rep"] == rep].copy()
            data_rep = data_rep[~data_rep["label_false_pos"]].copy()

            data_rep["start_peak"], data_rep["end_peak"]  = data_rep["peak_location_actual"].astype(int), data_rep["peak_location_actual"].astype(int)
            data_rep["end_peak"] = data_rep["end_peak"] +1
            data_rep["name"], data_rep["strand"] = """'.'""", """'.'"""
            data_rep["signal_value"] = data_rep["normalized_score_peak"].astype(int)
            data_rep["p_value"], data_rep["q_value"], data_rep["peak"] = 0,0,-1

            data_location_bed = data_rep[
                ["seqname", 'start_peak', 'end_peak', 
                "name", "normalized_score_peak", "strand", 
                "signal_value","p_value", "q_value","peak" ]
                ].copy()


            data_location_bed.sort_values(["seqname",  'start_peak', 'end_peak'], inplace = True)
            data_location_bed.copy().drop_duplicates(inplace= True)
            
            data_location_bed.to_csv(
                os.path.join(self.WP_directory, f"PeakLocations_{self.designfactor_value}_{rep}.narrowPeak"),
                sep  ="\t", 
                header = None, 
                index =False
                )


    def peak_ranges2narrowPeak(self, data):
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()

        for rep in data[f"{tagged_factor}_rep"].unique():
            data_rep = data[data[f"{tagged_factor}_rep"] == rep].copy()
            data_rep = data_rep[~data_rep["label_false_pos"]].copy()
            data_rep["peak_location_start"], data_rep["peak_location_stop"]  = data_rep["peak_location_start"].astype(int) , data_rep["peak_location_stop"].astype(int)


            data_rep["name"], data_rep["strand"] = """'.'""", """'.'"""
            data_rep["signal_value"] = data_rep["normalized_score_peak"].astype(int)
            data_rep["p_value"], data_rep["q_value"], data_rep["peak"] = 0,0,-1

            data_ranges_bed = data_rep[
                ["seqname", 'peak_location_start', 'peak_location_stop',
                "name", "normalized_score_peak", "strand", 
                "signal_value","p_value", "q_value","peak" ]].copy()
            data_ranges_bed.sort_values(["seqname", 'peak_location_start', 'peak_location_stop'], inplace = True)
            data_ranges_bed.copy().drop_duplicates(inplace=True)
            

            data_ranges_bed.to_csv(os.path.join(self.WP_directory, f"PeakSpans_{self.designfactor_value}_{rep}.narrowPeak"), sep  ="\t", header = None, index =False)
            
            
    def save_data(self, data, data_filter, peaks2gtf, keep_false_pos = False):
        designfactor_dict, tagged_factor, other_factor = self.get_tag_factor()
        filename_add = "_include_false_pos"
        if not keep_false_pos:
            data = data[~data["label_false_pos"]].copy()
            peaks2gtf = peaks2gtf[~peaks2gtf["label_false_pos"]].copy()
            filename_add = ""
        
        self.peak_locations2narrowPeak(data)
        self.peak_ranges2narrowPeak(data)
        for df, filename in zip([
            data, data_filter, peaks2gtf], ["all_taggedpeaks", "tagged2untagged_peaks", "peaks2annotation"]
                            ):
            df.to_csv(os.path.join(self.WP_directory,f"WP_out_{self.designfactor_value}_{filename}{filename_add}.csv"), index=False)
            
            
            
