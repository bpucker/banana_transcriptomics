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



### References ###


