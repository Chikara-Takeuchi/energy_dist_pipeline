import os
import sys

import gc
import pickle
import pandas as pd
import numpy as np
import scanpy as sc
from tqdm import tqdm
from collections import defaultdict # Import defaultdict

def load_annotation(annotation_file_path):
    """
    Load annotation file. check format(tab, or comma-delimed) and load with appropriate format.
    """
    print("Load csv formated annotation file")
    annotation_df = pd.read_csv(annotation_file_path,index_col=None)

    if ("guide_id" not in annotation_df.columns) | ("type" not in annotation_df.columns):
        print("required columns(guide_id,type) is not detected in the columns. Try loading as tab-delimed files")
        annotation_df = pd.read_csv(annotation_file_path,sep="\t",index_col=None)
        if ("guide_id" not in annotation_df.columns) | ("type" not in annotation_df.columns):
            print("[Error] Annotation file format. \n",
                  "Are required columns included in the annotation files? is this a comma or tab-delimed file?")
            sys.exit(1)
    print("annotation file loaded")
    print("--------")
    gRNA_type_list = [
                    "safe-targeting",
                    "non-targeting",
                    "targeting",
                    "positive control",
                    "negative control",
                    "variant"]
    for gRNA_type in gRNA_type_list:
        gRNA_type_anno = annotation_df[annotation_df["type"]==gRNA_type]
        print(f"{gRNA_type}: {gRNA_type_anno.shape[0]} gRNAs")
    return annotation_df

def check_annotation_gRNA_table(annotation_df,gRNA_dict):
    if not "guide_id" in annotation_df.columns:
        print("annotation df does not contains [guide_id] column name. Please check the format")
        sys.exit(1)
    else:
        annotation_gRNA_names = annotation_df["guide_id"].values
        gRNA_name = list(gRNA_dict.keys())
        overlap_gRNA = np.intersect1d(annotation_gRNA_names,gRNA_name)
        annotation_only = np.setdiff1d(annotation_gRNA_names,gRNA_name)
        gRNA_dict_only = np.setdiff1d(gRNA_name,annotation_gRNA_names)
        
        print(f"{len(overlap_gRNA)} gRNAs are found in both annotation and gRNA_dict")
        print(f"{len(annotation_only)} gRNAs are found only in annotation (possibly dropout)")
        print(f"{len(gRNA_dict_only)} gRNAs are found only in gRNA_dict")
        
        annotation_only_df = pd.DataFrame(columns=["gRNA_name","status"])
        annotation_only_df["gRNA_name"] = annotation_only
        annotation_only_df["status"] = "annotation_only"
        
        gRNA_dict_only_df = pd.DataFrame(columns=["gRNA_name","status"])
        gRNA_dict_only_df["gRNA_name"] = gRNA_dict_only
        gRNA_dict_only_df["status"] = "gRNA_dict_only"
        
        output_df = pd.concat([annotation_only_df,gRNA_dict_only_df])
        return output_df        
    
def get_gRNA_region_dict(annotation_df,gRNA_dict,
                         concatenate_key="intended_target_name",non_target_key="non-targeting"):
    """
    Creates a dictionary mapping target transcript names to a unique list of associated gRNA protospacer IDs.

    Args:
        annotation_df (pandas.DataFrame): DataFrame containing annotation information
        gRNA_dict (dict): A dictionary where keys are gRNA protospacer IDs. The values are not used in this function.

    Returns:
        dict: A dictionary where keys are target transcript names and values are NumPy arrays
              containing unique gRNA protospacer IDs associated with that transcript.
    """
    
    gRNA_region_dict_tmp = {}
    gRNA_region_dict_tmp[non_target_key] = []
    for index,row in annotation_df.iterrows():
        if row["guide_id"] in gRNA_dict.keys():
            if row["type"] == "non-targeting":
                gRNA_region_dict_tmp[non_target_key] += [row["guide_id"]]
            else:
                if row[concatenate_key] in gRNA_region_dict_tmp.keys():
                    gRNA_region_dict_tmp[row[concatenate_key]] += [row["guide_id"]]
                else:
                    gRNA_region_dict_tmp[row[concatenate_key]] = [row["guide_id"]]
    for key in gRNA_region_dict_tmp.keys():
        gRNA_region_dict_tmp[key] = np.unique(gRNA_region_dict_tmp[key])
    
    return gRNA_region_dict_tmp

def load_files(input_file, sgRNA_file, pca_file, dict_file, obsm_key="X_pca", overwrite=False):
    """
    Loads input and sgRNA files, returns PCA data and a gRNA dictionary.
    Reads from intermediate files if they exist, otherwise generates and saves them.
    If overwrite=True, forces regeneration and overwriting of intermediate files.

    Args:
        input_file (str): Input AnnData file (.h5ad).
        sgRNA_file (str): sgRNA count data (DataFrame in pickle format).
        pca_file (str): Path to save/load the PCA results pickle file.
        dict_file (str): Path to save/load the gRNA dictionary pickle file.
        obsm_key (str, optional): Key for PCA data within the AnnData object. Defaults to "X_pca".
        overwrite (bool, optional): Whether to force overwrite intermediate files. Defaults to False.

    Returns:
        tuple: (pd.DataFrame, dict): PCA data (X) and gRNA dictionary (gRNA_dict).
               Returns (None, None) or (X, None) if critical errors occur.
    """
    print("--- Processing PCA data ---")
    # Load from AnnData if pca_file doesn't exist or overwrite is True
    if not os.path.isfile(pca_file) or overwrite:
        if not os.path.isfile(pca_file):
             print(f"PCA file '{pca_file}' not found. Loading PCA from AnnData.")
        else:
             print(f"Overwrite specified. Loading PCA from AnnData (will overwrite '{pca_file}').")

        print(f"Loading input file '{input_file}'...")
        try:
            input_data = sc.read_h5ad(input_file)
            input_data_index = list(input_data.obs.index)
            print(f"Extracting PCA data from '{obsm_key}'...")
            # Copying is safer when creating the DataFrame
            X = pd.DataFrame(input_data.obsm[obsm_key].copy(), index=input_data_index)
            print(f"Saving PCA data to '{pca_file}'...")
            X.to_pickle(pca_file)
            print("Releasing AnnData object from memory.")
            del input_data # Release memory early
            gc.collect() # Garbage collection
        except FileNotFoundError:
            print(f"Error: Input file '{input_file}' not found.")
            return None, None
        except KeyError:
            print(f"Error: Key '{obsm_key}' not found in the AnnData object.")
            return None, None
        except Exception as e:
            print(f"An error occurred during AnnData processing or PCA saving: {e}")
            return None, None
    else:
        print(f"Loading existing PCA file '{pca_file}'.")
        try:
            X = pd.read_pickle(pca_file)
            input_data_index = list(X.index) # Get index for later use
        except FileNotFoundError:
            print(f"Error: Specified PCA file '{pca_file}' not found, although it was expected.")
            return None, None
        except Exception as e:
            print(f"An error occurred while loading PCA file '{pca_file}': {e}")
            return None, None


    print("\n--- Processing gRNA dictionary ---")
    # Create dictionary from sgRNA data if dict_file doesn't exist or overwrite is True
    if not os.path.isfile(dict_file) or overwrite:
        if not os.path.isfile(dict_file):
            print(f"Dictionary file '{dict_file}' not found. Creating dictionary from sgRNA data.")
        else:
            print(f"Overwrite specified. Creating dictionary from sgRNA data (will overwrite '{dict_file}').")

        print(f"Loading sgRNA file '{sgRNA_file}'...")
        try:
            sgRNA_data = pd.read_pickle(sgRNA_file)
        except FileNotFoundError:
            print(f"Error: sgRNA file '{sgRNA_file}' not found.")
            return X, None
        
        except MemoryError:
            print(f"MemoryError: Not enough memory to load sgRNA file '{sgRNA_file}'.")
            return X, None
        
        except Exception as e:
            print(f"An error occurred while loading sgRNA file '{sgRNA_file}': {e}")
            return X, None

        print("Transposing sgRNA data and filtering by cells present in input data...")
        try:
            sgRNA_data = sgRNA_data.T
            # Keep only cells present in the index of X
            sgRNA_data = sgRNA_data[sgRNA_data.index.isin(input_data_index)]
        except MemoryError:
            print(f"MemoryError: Not enough memory to load sgRNA file '{sgRNA_file}'.")
            return X, None
        
        except Exception as e:
            print(f"An error occurred during sgRNA data transpose or filtering: {e}")
            del sgRNA_data
            gc.collect()
            return X, None


        print("Extracting non-zero count cell-gRNA pairs...")
        # Using scipy.sparse would be more efficient if the matrix is sparse
        # Processing with numpy here, similar to the original code, but potential for improvement if data is sparse
        try:
             non_zero_rows, non_zero_cols = np.where(sgRNA_data.values != 0)
        except Exception as e:
            print(f"An error occurred while finding non-zero elements in sgRNA data: {e}")
            del sgRNA_data
            gc.collect()
            return X, None


        if len(non_zero_rows) == 0:
            print("Warning: No non-zero counts found in the sgRNA data. Returning an empty dictionary.")
            gRNA_dict = {}
        else:
            gRNA_name_list = sgRNA_data.columns.to_list()
            cell_name_list = sgRNA_data.index.to_list() # Use index after filtering

            print("Creating gRNA dictionary (improved method)...")
            # Build the dictionary efficiently using defaultdict
            gRNA_dict_default = defaultdict(list)
            num_pairs = len(non_zero_rows)
            try:
                for i in tqdm(range(num_pairs), desc="Processing gRNA pairs"):
                    row_idx = non_zero_rows[i]
                    col_idx = non_zero_cols[i]
                    cell_name = cell_name_list[row_idx]
                    gRNA_name = gRNA_name_list[col_idx]
                    gRNA_dict_default[gRNA_name].append(cell_name)

                # Convert defaultdict back to a regular dict
                gRNA_dict = dict(gRNA_dict_default)
                print(f"gRNA dictionary creation complete. Found {len(gRNA_dict)} types of gRNAs.")
            except Exception as e:
                print(f"An error occurred during gRNA dictionary creation loop: {e}")
                # Clean up intermediate variables even if loop fails
                del sgRNA_data, non_zero_rows, non_zero_cols
                if 'gRNA_dict_default' in locals(): del gRNA_dict_default
                gc.collect()
                return X, None # Return X, but indicate gRNA dict failed


        print(f"Saving gRNA dictionary to '{dict_file}'...")
        try:
            with open(dict_file, mode='wb') as fo:
                pickle.dump(gRNA_dict, fo)
        except IOError as e:
            print(f"Error: Failed to write gRNA dictionary file '{dict_file}': {e}")
        except Exception as e:
            print(f"An unexpected error occurred while saving the gRNA dictionary: {e}")


        print("Releasing sgRNA data from memory.")
        del sgRNA_data # Release memory
        del non_zero_rows, non_zero_cols # Delete variables that are no longer needed
        if 'gRNA_dict_default' in locals():
             del gRNA_dict_default
        gc.collect() # Garbage collection

    else:
        print(f"Loading existing gRNA dictionary file '{dict_file}'.")
        try:
            with open(dict_file, mode='br') as fi:
                gRNA_dict = pickle.load(fi)
            print(f"gRNA dictionary loaded successfully. Found {len(gRNA_dict)} types of gRNAs.")
        except FileNotFoundError:
            print(f"Error: Specified gRNA dictionary file '{dict_file}' not found, although it was expected.")
            return X, None # Return X, but indicate gRNA dict failed
        except Exception as e:
            print(f"Error: An error occurred while loading gRNA dictionary file '{dict_file}': {e}")
            return X, None

    print("\n--- Processing finished ---")
    # Removed potentially unnecessary gc.collect() here
    return X, gRNA_dict



def get_unique_list(seq):
    return_list = []
    seen = []
    for x in seq:
        if x in seen:
            continue
        elif (x[1],x[0]) in seen:
            continue
        else:
            seen.append(x)
    return seen

def extract_gene_name(entry):
    if entry.startswith('OR'):
        return entry.split("-")[0]
    else:
        return entry.split("_")[0]
    
def extract_transcript_name(entry):
    gene_name = ""
    transcript_name = ""
    if entry.startswith('OR'):
        gene_name = entry.split("-")[0]
    else:
        gene_name = entry.split('_')[0]
        transcript_name = entry.split('-')[-2:][0]
    if len(entry.split('_')) <3:
        transcript_name = gene_name
    else:
        transcript_name = gene_name + ':' + transcript_name
    return transcript_name