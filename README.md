# Banana Transcriptomics #

These scripts are associated with banana transcriptomic studies.

### Co-expression analysis ###

This script is based on a previously developed script [MtMYB](https://github.com/bpucker/MtMYBs). Pairwise co-expression is analyses for a given list of gene IDs against all gene IDs in the expression data file.

```
python coex_analysis.py --candidates <FILE> --out <DIR> --exp <FILE> --anno <FILE>

Mandatory:
  --candidates    STR         Sequence ID input file
  --out           STR         Output folder
  --exp           STR         Expression data file
  --anno          STR         Annotation file
```

`--candidates` specifies a text file that contains gene IDs of interest. One ID per line is expected.

`--out` specifies a output folder for the result files. This folder will be created if it does not exist already.

`--exp` specifies an expression data file. The first row contains the sample IDs and the first column contains the gene IDs.

`--anno` specifies an annotation text file with gene IDs in the first column and an annotation text in the second column.


### Construct Cytoscape input file ###

This script generates an input file for cytoscape to visualize the interaction of transcription factors and target genes. We investiagted numerous transcription factor families including MYB, bHLH, WD40, WRKY, MADS-box, ARF, bZIP, NAC, and Hsf and their target genes in the flavonoid biosynthesis. However, this script could be used to investigate the regulation of other pathways by these transcription factor families.

```
python construct_cytoscape_input.py --genes <FILE> --out <FILE> --coexp <DIR>

Mandatory:
  --genes    STR         Sequence ID input file
  --out      STR         Output file
  --coexp    STR         Coexpression result folder
  
  optional:
  --tf       STR         Transcription factor file
  --r        FLOAT       Correlation coefficient cutoff
  --p        FLOAT       p-value cutoff
```

`--genes` specifies a text file that contains gene IDs of interest. One ID per line is expected.

`--out` specifies a text file that contains the output. This file can be imported into cytoscape for visualization. Columns: Target genes, transcription factor, correlation (r), and transcription factor class.

`--coexp` specifies a folder which contains co-expression analysis result files. Data of all files will be loaded, but only the genes specified in the genes file and the TF file will be considered for the output dataset.

`--tf` specifies a text file which describes the transcription factors. The first column contains the gene ID, while the second column contains the name which also indicates the family.

`--r` specifies a correlation coefficient cutoff value. Only gene pairs with higher correlation than this cutoff will be considered for the output file.

`--p` specifies a p-value of the correlation coefficient of a gene pair. Only gene pairs with lower p-values than this cutoff will be considered for the output file.



### References ###


