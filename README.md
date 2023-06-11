# cell sub-crosstalk pair
Defining and identifying cell sub-crosstalk pairs for characterizing cell-cell communication patterns
![image](https://github.com/chenxing-zhang/cell-sub-crosstalk-pair/blob/main/schematic_diagram.png)

## Descriptions
We define a cell sub-crosstalk pairs of two cell types contains two cell subgroups with strong and similar crosstalk signals and identify cell sub-crosstalk pairs based on coupled non-negative matrix factorization and coupled consensus clustering.

## Preparations
**1. Datasets**  
Demo datasets are saved in `"..//Demo_data"`, including scRNA-seq matrix of mouse endothelial cells `"mat_F002_Endo"` and microglia cells `"mat_F002_Micro"` and known ligand-receptor pairs (`"PairsLigRec_mouse.txt"` and `"PairsLigRec_human.txt"` ).

**2. Package requirements**  
Python packages `sklearn`, `numpy`, `pandas`, `statsmodels` are required

## Codes
The scipts of identifying cell sub-crosstalk pairs are in `“..//Code//identify_CSCPs.py”`, which can be executed using following commands:
```
python identify_CSCPs.py [-h] [--lrPath LRPATH] [--lrName LRNAME]
                         [--matPath MATPATH] [--matNameL MATNAMEL] [--matNameR MATNAMER]
                         [--K K] --outPath OUTPATH
Arguments:
  -h, --help           show this help message and exit
  --lrPath LRPATH      the path of Ligand-Receptor dataset
  --lrName LRNAME      the name of Ligand-Receptor dataset
  --matPath MATPATH    the path of expression matrix
  --matNameL MATNAMEL  the name of expression matrix(XL) of sender cell type
  --matNameR MATNAMER  the name of expression matrix(XR) of sender cell type
  --K K                the number of CSCPs. default=2
  --outPath OUTPATH    the path to save result
```

Identifying cell sub-crosstalk pairs between mouse endothelial cells and microglia cells:
```
python identify_CSCPs.py --lrPath ..//Demo_data --lrName PairsLigRec_mouse.txt --matPath ..//Demo_data --matNameL mat_F002_Endo --matNameR mat_F002_Micro --outPath ..//Demo_data
```

## Reference
Zhang C, Hu Y, Gao L. Defining and identifying cell sub-crosstalk pairs for characterizing cell-cell communication patterns. (under review)