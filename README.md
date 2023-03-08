# cell sub-crosstalk pair
Defining and identifying cell sub-crosstalk pairs for characterizing cell-cell communication patterns
![image](https://github.com/chenxing-zhang/cell-sub-crosstalk-pair/blob/main/schematic_diagram.png)

## Descriptions
We define a cell sub-crosstalk pairs of two cell types contains two cell subgroups with strong and similar crosstalk signals and identify cell sub-crosstalk pairs based on coupled non-negative matrix factorization and coupled consensus clustering.

## Preparations
**1. Datasets**  
Demo datasets are saved in `".//Demo_data"`, including scRNA-seq matrix of mouse endothelial cells `(“.//Demo_data//mat_F002_Endo”)` and microglia cells `(“.//Demo_data//mat_F002_Micro”)` and Known ligand-receptor pairs `(“.//Demo_data//PairsLigRec_simple_mouse.txt”)`.

**2. Package requirements**  
Python packages `sklearn`, `numpy`, `pandas`, `statsmodels` are required

## Codes
The code for identifying cell sub-crosstalk pairs is in `(“.//Code//identifying_cell_sub-crosstalk_pairs.py”)`, which includes the the following steps:
1. Reading ligand-receptor pairs and single cell RNA-seq data (function `read_csv`).
2. Preprocessing single cell RNA-seq data (function `preprocessing_expression_matrix`).
3. Constructing ligand-receptor matrix `A` (function `construct_ligand_receptor_matrix`).
4. Identifying candidates of cell sub-crosstalk pairs by coupled non-negative matrix factorization (function `identifying_cellSubCrossTalkPairs_candidates`).
5. Constructing coupled concensus matrix (function `construct_coupled_concensus_matrix`). 
6. Merging candidates by hierarchical clustering (function `hierarchical_clustering`).
7. Adding label of cell sub-crosstalk pairs (function `create_cellSubCrossTalkPairs_label`).

## Reference
Zhang C, Hu Y, Gao L. Defining and identifying cell sub-crosstalk pairs for characterizing cell-cell communication patterns. (under review)