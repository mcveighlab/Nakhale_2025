# Nakhale_2025
Analysis code for bioinformatics datasets reported in Nakhale et al. 2025

This code was used to analyse small RNA sequencing reads from NGS mice carrying experimental _Dirofilaria immitis_ infections. 
Plasma was recovered from infected mice and uninfected (sham injected) controls at days 5, 28 and 70 following infection.
RNA was extracted from plasma samples, and sequenced by Illumina small RNA methods.
Sequencing reads were analysed for the presence of micro (mi)RNAs using the miRDeep2_subtractive.sh code, generating counts for host and parasite derived miRNAs.
Host miRNA counts were then analysed using the DESeq2.R code to identify statistically significant differentially expressed miRNAs, as reported in the paper.
