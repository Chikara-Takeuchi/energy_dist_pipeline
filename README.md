
# TF Perturb-Seq Energy Distance Pipeline

The TF Perturb-Seq energy distance pipeline is designed to perform the following tasks:

1.  Filter outlier gRNAs from Perturb-Seq data.
2.  Identify significant perturbations (hits) from Perturb-Seq screens using energy distance and permutation tests.
3.  Identify perturbations with similar effects to characterize distinct perturbation phenotypes.

Inputs are:

* An `.h5ad` file containing processed single-cell data, including a PCA embedding (typically found in `.obsm['X_pca']`).
* A cell-by-gRNA count matrix (saved in Python pickle format).
* (Optional) gRNA information table(csv format) column should include `guide_id`,`intended_target_name`,`type`,`spacer` defined in the [IGVF DACC](https://github.com/IGVF-DACC/checkfiles/blob/dev/src/schemas/table_schemas/guide_rna_sequences.json)

gRNA information table should located in `OUTPUT_FOLDER`/`annotation_file` specified in the config file.



## Run TF Perturb-Seq Energy Distance Pipeline

**First, create the necessary configuration files. See the "Configuration Files" section below for details.**

### With Container (Recommended)

1.  Modify the `run_container_step0.sh` and `run_container_step1_2.sh` file as follows:
    ```bash
    CONTAINER_PATH=[PATH TO THE .sif FILE FOR THE ENERGY DISTANCE CONTAINER]
    CONFIG_PATH=[PATH TO YOUR config.json]
    BIN_PATH=[PATH TO THE bin FOLDER OF THIS PIPELINE]
    ```

2.  (Optional) Format gRNA information. if you already have gRNA information table with IGVF DACC format, skip this step.
    ```bash
    sbatch run_container_step0.sh
    ```

3. Then, run the pipeline using:
    ```bash
    sbatch run_container_step1_2.sh
    ```

### Without Container

1.  Create a Python environment that satisfies the dependencies listed in `requirements.txt`.
2.  Modify the `run_step0.sh` and `run_step1_2.sh` file as follows:
    ```bash
    CONFIG_PATH=[PATH TO YOUR config.json]
    BIN_PATH=[PATH TO THE bin FOLDER OF THIS PIPELINE]

    # Load your environment (example for Conda)
    source activate [YOUR_ENVIRONMENT_NAME]
    ```
3.  Run the pipeline using:
    ```bash
    sbatch run_step0.sh
    sbatch run_step1_2.sh
    sbatch run_step3.sh
    ```

## Descriptions of each script

This pipeline consists of the following Python scripts:

### `0_preprocess.py`

* **Purpose**: Performs initial setup for the pipeline and prepares the annotation file.
* **Function**:
    * Reads the reference gRNA annotation file (e.g., `Hon_sgRNA_index_dacc_annot_reference.csv`), along with lists of negative controls and non-targeting controls.
    * For each gRNA, identifies its target gene/transcript and determines its source type (target, non-targeting, negative control, positive control).
    * Outputs a formatted annotation file (named according to `annotation_file` in `config.json`) for use by subsequent scripts.

### `1_filtereing_gRNA.py`

* **Purpose**: Performs quality filtering of gRNAs. Specifically, it assesses the consistency of effects among multiple gRNAs targeting the same region, filtering out outlier gRNAs and non-targeting gRNAs that may not represent a true baseline.
* **Function**:
    * **Targeting gRNA Filtering**:
        * Calculates Energy Distance between pairs (or subsets of pairs) of gRNAs targeting the same region.
        * Runs the DISCO (Distance Components) test to evaluate if the cell state distributions (in PCA space) induced by different gRNAs for the same target are statistically different.
        * Based on the ranks from the intra-region Energy Distance calculations, uses a hypergeometric test to compute a p-value for each gRNA, identifying potential outliers.
    * **Non-targeting gRNA Filtering**:
        * Calculates pairwise Energy Distances among non-targeting gRNAs.
        * Uses K-means clustering on the resulting distance matrix to identify and filter out non-targeting gRNAs that do not belong to the majority cluster, considering them outliers.
    * Outputs the filtering results for targeting and non-targeting gRNAs into separate files (`targeting_outlier_table`, `non_targeting_outlier_table`).

### `2_e_distance_nontargeting.py`

* **Purpose**: Evaluates the strength and statistical significance of the perturbation effect for each target region (using filtered gRNAs).
* **Function**:
    * Identifies the cell populations corresponding to each target region, using only the non-outlier gRNAs identified in Step 1.
    * Calculates the Energy Distance between the cell population for each target region and multiple, randomly sampled sets of non-targeting control cells (using filtered non-targeting gRNAs).
    * Performs a permutation test for each comparison to assess the statistical significance (p-value) of the observed Energy Distance.
    * Calculates the mean Energy Distance and mean p-value across the multiple non-targeting background comparisons for each target region and outputs the results to a table (`edist_pvalue_table`).

### `2_1_Plot_figure.py`

* **Purpose**: Generates visualizations of intermediate and final results for quality control and parameter selection assessment.
* **Function**:
    * Plots histograms showing the count of gRNAs per target region categorized as dropped, outlier, or retained after the filtering in Step 1.
    * Plots the distribution of Energy Distance versus p-value calculated in Step 2 (scatter plot and density plot).
    * Generates heatmaps showing the number of target regions that pass significance thresholds under varying combinations of p-value and Energy Distance cutoffs, aiding in the selection of appropriate cutoffs for Step 3.

### `3_e_distance_among_regions.py`

* **Purpose**: Compares target regions that showed significant perturbation effects in Step 2, clusters regions with similar effects, and generates an embedding for visualization.
* **Function**:
    * Selects target regions with significant perturbation effects based on the results from Step 2 (`edist_pvalue_table`) and the cutoffs specified in `config_clustering.json` (p-value and Energy Distance).
    * Calculates a pairwise Energy Distance matrix between all selected significant target regions.
    * Performs clustering (currently primarily Affinity Propagation, using the negative of the distance matrix) to group regions with similar perturbation profiles.
    * Generates a 2D t-SNE embedding based on the pairwise distance matrix.
    * Outputs the pairwise distance matrix (`edist_target_by_target_matrix`) and the t-SNE embedding coordinates with cluster assignments (`edist_embedding_info`).


## Configuration Files

This pipeline uses two configuration files in JSON format: `config.json` for the main pipeline steps and `config_clustering.json` specifically for the final clustering step.

### `config.json`

This file contains parameters controlling input/output, gRNA filtering, and permutation testing.

#### `output_file_name_list`

This section defines the paths and filenames for various output and intermediate files generated throughout the pipeline. This part is not ususally needed to be edited.

* **`OUTPUT_FOLDER`**: (String) The main directory where all output files and subdirectories (like figures) will be saved. 

* **`targeting_outlier_table`**: (String) Filename for the table storing outlier p-values for *targeting* gRNAs. Generated by `1_filtereing_gRNA.py` based on DISCO and hypergeometric tests. 

* **`non_targeting_outlier_table`**: (String) Filename for the table storing outlier status for *non-targeting* gRNAs. Generated by `1_filtereing_gRNA.py` based on pairwise energy distance and K-means clustering of non-targeting gRNAs. 

* **`edist_pvalue_table`**: (String) Filename for the table storing the results of the energy distance calculation and permutation testing between each target region and the non-targeting background. Generated by `2_e_distance_nontargeting.py`.

* **`edist_target_by_target_matrix`**: (String) Filename for the matrix containing pairwise energy distances between all *significant* target regions (selected based on cutoffs in `config_clustering.json`). Generated by `3_e_distance_among_regions.py`. This matrix is used as input for clustering and t-SNE embedding.

* **`edist_embedding_info`**: (String) Filename for the table storing the t-SNE embedding coordinates (x, y) and cluster assignments for the significant target regions. Generated by `3_e_distance_among_regions.py`.

* **`pca_table`**: (String) Filename for the intermediate pickle file storing the PCA result (cell embeddings) extracted from the input AnnData file. This is handled by the `util_functions.load_files` function called in `1_filtereing_gRNA.py`, `2_e_distance_nontargeting.py`, and `3_e_distance_among_regions.py`. If the file exists, it's loaded; otherwise, it's generated from the `h5ad_file`.

* **`gRNA_dict`**: (String) Filename for the intermediate pickle file storing a dictionary mapping gRNA names to the list of cell barcodes associated with each gRNA. This is handled by the `util_functions.load_files` function called in `1_filtereing_gRNA.py`, `2_e_distance_nontargeting.py`, and `3_e_distance_among_regions.py`. If the file exists, it's loaded; otherwise, it's generated from the `sgRNA_file`.

#### `input_data`
This section defines the paths to the primary input data files.

* **`annotation_file`**: 
    configs related to annotation file.
    * **`file_path`**
    (String) Path to the processed gRNA annotation csv table. Column name should include `guide_id`,`intended_target_name`,`type`,`spacer`. additional columns are acceptable.
    * **`concatenate_key`**
    (String) the columns name of the annotation file to concatenate gRNAs. It is usually the region name (such as "intended_target_name").

* **`h5ad_file`**: 
    configs related to input h5ad file.
    * **`file_path`**
    (String) Path to the input AnnData file (`.h5ad`). This file should contain the gene expression data (though not directly used by these scripts) and, importantly, the precomputed PCA results stored in `.obsm['X_pca']`. Used by `util_functions.load_files` (called in scripts 1, 2, 3).
     * **`obsm_key`**
     (String) key in `obsm` of h5ad file for the e-dist calculation.

* **`sgRNA_file`**: 
    * **`file_path`**
    (String) Path to the sgRNA count matrix stored in pickle format (`.pkl`). This file should contain a DataFrame where rows are cells and columns are gRNAs (or vice-versa, the script transposes it if necessary). Used by `util_functions.load_files` (called in scripts 1, 2, 3).

#### `gRNA_filtering`

Parameters controlling the gRNA filtering process in `1_filtereing_gRNA.py`.
* **`perform_targeting_filtering`**
(bool) Whether this pipeline filter targeting gRNAs or not.
* **`perform_nontargeting_filtering`**
(bool) Whether this pipeline filter non-targeting gRNAs or not.
* **`threshold_gRNA_num`**
(Integer) If the number of gRNAs per region is less than this threshold, all of the combinations are used to find outlier gRNAs. if the number of gRNAs is more than this threshold, combinations within `combi_count` gRNAs are used.
* **`combi_count`**: (Integer) Used when calculating pairwise energy distances *within* a target region's gRNAs. If a target region has more than 6 associated gRNAs, this parameter defines the size of the subsets of gRNAs to compare (e.g., compare a subset of size `k` vs. another subset of size `combi_count - k`). If the region has 6 or fewer gRNAs, all possible pairs of disjoint non-empty subsets are compared. 

* **`total_permute_disco`**: (Integer) The number of permutations to perform for the DISCO (Distance Components) test. The DISCO test assesses if the distributions of cells associated with different gRNAs for the *same* target region are significantly different from each other.

* **`batch_num_basic`**: (Integer) The base batch number used in the DISCO test calculation. The actual batch size might be adjusted downwards based on the total number of cells involved to manage memory, particularly for GPU calculations.

#### `permutation_test`

Parameters controlling the permutation testing in `2_e_distance_nontargeting.py`, where each target region is compared against the non-targeting background.

* **`permute_per_bg`**: (Integer) The number of permutations performed for *each* comparison between a target region and one instance of the randomly sampled non-targeting background. 

* **`num_of_bg`**: (Integer) The number of independent non-targeting background cell sets to generate. Each target region is compared against each of these background sets, and the resulting p-values and distances are averaged to get a more robust estimate.

* **`non_target_pick`**: (Integer) The number of cells to randomly sample (from the pool of all cells associated with *filtered* non-targeting gRNAs) to create each background set specified by `num_of_bg`. Therefore, total number of the non-targeting cells used for background is `num_of_bg` * `non_target_pick`.

* **`batch_num_basic`**: (Integer) The base batch number used for the energy distance permutation test calculations (target vs. non-targeting). Similar to the DISCO test, the actual batch size might be adjusted based on the number of cells in the target region to manage computational resources.

#### `aggregate`

Parameters controlling the aggregation and final comparison step in `3_e_distance_among_regions.py`.

* **`downsampling_maximum`**: (Integer) The maximum number of cells per target region to use when calculating the pairwise energy distance matrix between *significant* target regions. If a target region (after filtering outlier gRNAs) has more cells than this value, its cells are randomly downsampled to this number before calculating the distance to other target regions. This helps to standardize comparisons and manage computational load.

## `config_clustering.json`

This file contains parameters specifically for selecting significant perturbations and clustering them in `3_e_distance_among_regions.py`.

#### `cutoff`

Defines the thresholds applied to the results from `2_e_distance_nontargeting.py` (stored in `edist_pvalue_table`) to select target regions considered to have a significant perturbation effect for further analysis (pairwise comparison, clustering, embedding).

* **`pval_cutoff`**: (Float) The maximum *mean* p-value (averaged across the `num_of_bg` comparisons) allowed for a target region to be considered significant. Regions with a p-value *below* this threshold are selected.
* **`distance_cutoff`**: (Float) The minimum *mean* energy distance (averaged across the `num_of_bg` comparisons) required for a target region to be considered significant. Regions with a distance *above* this threshold are selected. A region must satisfy both `pval_cutoff` and `distance_cutoff` to be included in the final pairwise comparison and clustering.

#### `clustering`

Parameters defining the clustering method for grouping significant target regions based on their pairwise energy distance matrix.

* **`method`**: (String) ["Affinity"]
Specifies the clustering algorithm. Although this parameter exists, the current implementation in `3_e_distance_among_regions.py` appears hardcoded to use Affinity Propagation (`sklearn.cluster.AffinityPropagation`) with precomputed distances (the negative energy distance matrix) regardless of the value set here.


# Contact:
**Chikara Takeuchi:** \
chikara.takeuchi [at] utsouthwestern.edu \
\
**Gary Hon:** \
gary.hon [at] utsouthwestern.edu