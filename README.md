## Introduction
This repository contains analysis for inference of non-interaction data using [Intact](https://www.ebi.ac.uk/intact/interactomes). Using this non-interaction data, high confidence non-interaction bait-prey protein pairs are selected.

The [Snakemake](https://snakemake.github.io/) workflow is as follows:

![rules](docs/figures/rules.png)


## Output:
Expected output files are listed is `SnakeFile.smk` and are as follows:
- Scored interactions
  - `"work_folder/inferred_search_space/analysis/bias_reduced_ppis/threshold_negatome.csv"`: High confidence non-interactions  
  - `"work_folder/inferred_search_space/analysis/bias_reduced_ppis/cell_line/high_confidence.csv"`: High confidence interactions
- Cell line:
  - `"work_folder/plots/cell_line_prey.png"`: OR of prey being identified per cell line (Fisher's exact test) 
  - `"work_folder/inferred_search_space/analysis/bias_reduced_ppis/cell_line/threshold_negatome.csv"`: Low/high cell line prey detectability among Negatome   
  - `"work_folder/inferred_search_space/analysis/bias_reduced_ppis/cell_line/high_confidence.csv"`: Low/high cell line prey detectability among HCIs
  - `"work_folder/inferred_search_space/analysis/cell_line/bait_prior.csv"`: Marginalised prey probability given cell line.
- Subcellular localisation:
  - `"work_folder/plots/localisation_OR_y2h_ms.png",`: Plot of detectability of protein pairs based on mass spectrometry and yeast to hybrid methods 

## Config
Parameters are set in `config.yaml` as:
- `formated_ppi:` Path to formated cell line PPI file
- `cell_line_ppis:` Path to cell line annotated PPI file 
- `remove_single_publications:` Boolean if single PPI studies should be filtered out
- `min_total_tests:` Minimum tests per prey-cell line Fisher's exact test
- `min_total_observed:` Minimum observed tests per prey-cell line Fisher's exact test
- `part_or_difference_cutoff:` OR cutoff
- `pseudo_n:` Prior strength for PPI-probability 
- `HCL_frac:` Probability threshold defining HCI
- `id_pattern:` set to `"gene_name"` or `"uniprot_id"` 
- `localisation_file:` path to localisation data 
- `selected_cell_lines:` 
- `ms: ` list of mass spectrometry detection methods 
- `y2h: ` list of yeast 2 hybrid methods

## Setup
Download the *Intact* miTab and put in `data/intact`. Set `cell_line_ppis` to the output `CL_annotated_bait_prey.csv` from [cell line curated](https://github.com/JoelAAs/Cell_line_curated_PPI).

## Run
```commandline
snakemake -s SnakeFile.smk -c 5 --configfile config.yaml 
```