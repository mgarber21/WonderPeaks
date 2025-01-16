import os
import pandas as pd

def initiate_user_inputs(directory):
    """
    Reads and processes user inputs from a CSV file to create a dictionary of parameters.

    Parameters:
    - directory (str): Path to the directory containing the "NGS_user_inputs.csv" file.

    Returns:
    - dict: A dictionary containing processed user inputs.
    """
    # Read the "NGS_user_inputs.csv" file into a DataFrame and convert it to a dictionary
    df_inputs = pd.read_csv(os.path.join(directory, "NGS_user_inputs.csv"), index_col=0)["User Inputs"]
    user_inputs_dict = df_inputs.to_dict()

    # Parse the 'reads' input into a list, stripping whitespace from each entry
    user_inputs_dict["reads"] = [read.strip() for read in user_inputs_dict["reads"].split(",")]

    # Convert specific keys to None or ensure their data types are consistent
    for key in ["adapter_sequence", "include_transcripts_above", "include_transcripts_below"]:
        if user_inputs_dict[key] == "None":
            user_inputs_dict[key] = None
        elif isinstance(user_inputs_dict[key], int):
            user_inputs_dict[key] = str(user_inputs_dict[key])  # Convert integers to strings for consistency

    # Construct the full path to the genome fasta file
    if user_inputs_dict["genome_directory"] in user_inputs_dict["genome_fasta_file"]:
        user_inputs_dict["genome_fasta_path"] = user_inputs_dict["genome_fasta_file"]
    else:
        user_inputs_dict["genome_fasta_path"] = os.path.join(
            user_inputs_dict["genome_directory"], user_inputs_dict["genome_fasta_file"]
        )

    # Construct the full path to the genome annotations file
    if user_inputs_dict["genome_directory"] in user_inputs_dict["genome_annotations_file"]:
        user_inputs_dict["genome_annotations_path"] = user_inputs_dict["genome_annotations_file"]
    else:
        user_inputs_dict["genome_annotations_path"] = os.path.join(
            user_inputs_dict["genome_directory"], user_inputs_dict["genome_annotations_file"]
        )

    return user_inputs_dict

