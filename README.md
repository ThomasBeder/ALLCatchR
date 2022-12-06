# ALLCatchR

_The gene expression classifier ALLCatchR identifies B-precursor ALL subtypes and underlying developmental trajectories_

_ALLCatchR_ was developed to predict:
- 21 BCP-ALL molecular subtypes (BCL2/MYC, CDX2/UBTF, CEBP, DUX4, ETV6::RUNX1, ETV6::RUNX1-like, HLF, Hyperdiploid, iAMP21, IKZF1 N159Y, KMT2A, Low hypodiploid, MEF2D, Near haploid, NUTM1 PAX5 P80R, PAX5alt, Ph-like, Ph-pos, TCF3::PBX1, ZNF384)
- Associations to B lymphopoiesis stages based on gene set enrichment analyses 
- Blast Count percentage
- Immunophenotype
- Patient's sex

![image](ALLCatchR_workflow.png)

## Publication
Beder et al. 2022 [The gene expression classifier ALLCatchR identifies B-precursor ALL subtypes and underlying developmental trajectories across age groups]

## Installation
open RStudio
install devtools and follow the installion guide https://github.com/r-lib/devtools
```
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")
```
install ALLCatchR 
```
devtools::install_github("ThomasBeder/ALLCatchR")
```

## Quickstart
If Counts.file is left ```NA``` ten test samples are predicted
```
library(ALLCatchR)
allcatch(Counts.file = NA, ID_class = "symbol",sep = "\t")
```

## Run ALLCatchR
As input ALLCatchR requires a single text file in which the first column represent the genes and the other columns the count data for each sample
```
library(ALLCatchR)
allcatch(Counts.file = NA, ID_class = "symbol", sep = "\t")
# Counts.file: /path/to/your/count/data, if left empty a test
# ID_class: gene names can be either "symbol", "ensemble_ID" or	"entrez_ID"
# sep: seperator of the text file usually "\t", "," or ";"
```

## output
ALLCatchR writes a ```predictions.tsv``` file to your current working directory with the following columns:
- sample: Sample ID
- Score: BCP-ALL subtype prediction score
- Prediction: Predicted subtype
- Confidence: Confidence of subtype predictions, **IMPORTANT** predictions of samples labeled as unclassified here should be considered unclassified
- BlastCounts: Predicted Blast Count percentage
- Sex: Patient's sex predition
- Score_sex: Sex prediciton probability
- Immuno: Immunophenotype prediction
- ScoreImmuno: Immunophenotype prediciton probability
- 10-16: Enrichment analysis results (singscore, https://github.com/DavisLaboratory/singscore) to gene sets defined for B-cell progenitors
- 37-37: SVM linear predictions for individual subtypes
- 38-59: Results from gene set based nearest neighbor analysis
