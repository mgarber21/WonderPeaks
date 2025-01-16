from Helper_Functions import *
from itertools import combinations, permutations, groupby


def get_tag_factor(designfactors, metadata):
    designfactors_col = "_".join(designfactors)

    designfactor_dict =  dict()
    for designfactor in designfactors:
        factor_samples = metadata[designfactor].unique()
        if [factor_sample for factor_sample in factor_samples if "tag" in factor_sample.lower()]:
            # tagged factor is designfactor column to use for lowest level
            tagged_factor = designfactor
        else:
            other_factor = designfactor
        designfactor_dict[designfactor] = list(factor_samples)
    return designfactor_dict, tagged_factor, other_factor

def make_tagged_factor_rep(df, designfactors, metadata):
    designfactor_dict, tagged_factor, other_factor = get_tag_factor(designfactors, metadata)
    df[f"{tagged_factor}_rep"] = df.apply(lambda row: row[tagged_factor]+"_"+str(row["replicate"]),
                                    axis= 1)
    return df




def get_sample_combos(df, designfactors, metadata):
    
    designfactor_dict, tagged_factor, other_factor = get_tag_factor(designfactors, metadata)
    ## create dictionaries using sample metadata
    #  Note control will always be first in the dictionart because dictionaries are sorted alphabetically (control, test)
    dict_sample_types = {}
    dict_sample_types["control"] = [
        i for i in df[tagged_factor].unique() if "control" in i.lower() or "ctrl" in i.lower()
        ][0]
    dict_sample_types["test"]= [
        i for i in df[tagged_factor].unique() if not ("control" in i.lower() or "ctrl" in i.lower())
        ][0]
    
    dict_sample_types_reps = {}
    dict_sample_types_reps["control"] = [i for i in df[f"{tagged_factor}_rep"].unique()if "control" in i.lower() or "ctrl" in i.lower() or "untagged" in i.lower() or "input" in i.lower()]
    dict_sample_types_reps["test"] = [i for i in df[f"{tagged_factor}_rep"].unique() if not ("control" in i.lower() or "ctrl" in i.lower() or "untagged" in i.lower() or "input" in i.lower())]
    
    
    #  Note control will always be first because dictionaries are sorted
    

    
    #  Note control will always be first because dictionaries are sorted
    Sample_type_combos = list(combinations(df[f"{tagged_factor}_rep"].unique(), 2))
    
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


def other_factor4file_name(other_factor, df):
    other_factor_value = df[other_factor].unique()
    return other_factor_value


def metadata_connect(df, metadata, designfactors):
    metadata["replicate"] = metadata.groupby("_".join(designfactors)).cumcount() + 1
    df  = pd.merge(df, metadata, on = "basename")
    
    return df