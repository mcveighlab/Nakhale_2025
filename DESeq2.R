# set working directory
setwd("<path_to_working_directory>") 

# Install DESeq2 if required
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
} 

if (!requireNamespace("DESeq2", quietly = TRUE)) { 
  BiocManager::install("DESeq2") 
} 

# load DESeq2 into environment
library("DESeq2") 

# Load in count data generated during mirdeep2 analysis, adjust column labels
countData <- as.matrix(read.csv("<filename.csv>", row.names = "miRNA_name")) 
colnames(countData)<- sub("^X", "", colnames(countData))  # Remove "X" prefix from column names 

# Load in condition data. This is a two column text file:
# Left column titled SampleID, containing a descriptor of each individual sample
# Right column titled Condition, containing a descriptor of an individual condition e.g. timepoint, treated/control, age etc.
PHENO_DATA <- "PhenoData.csv" 
colData <- read.csv(PHENO_DATA, sep = ",", row.names = 1) 
colData$Condition <- as.factor(colData$Condition) 

# Create and run DESeq2 dataset 
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Condition) 
dds <- DESeq(dds) 

# Extract normalized counts 
normalized_counts <- counts(dds, normalized = TRUE) 

# Save normalized counts to a CSV file 
write.csv(normalized_counts, file = "Normalized_Counts.csv") 

# Define the specific comparisons you want to make 
specified_comparisons <- list( 
  c("D5_Inf", "D5_Sham"), 
  c("D28_Inf", "D28_Sham"), 
  c("D70_Inf", "D70_Sham"), 
  c("D5_Inf", "D28_Inf"), 
  c("D28_Inf", "D70_Inf"), 
  c("D5_Inf", "D70_Inf"), 
  c("D5_Sham", "D28_Sham"), 
  c("D28_Sham", "D70_Sham"), 
  c("D5_Sham", "D70_Sham") 
) 

# Perform specified comparisons 
results_list <- list() 
for (comp in specified_comparisons) { 
  contrast <- c("Condition", comp[1], comp[2]) 
  res <- results(dds, contrast = contrast) 

  # Save results to a file and list 
  result_name <- paste0(comp[1], "_vs_", comp[2], ".csv") 
  write.csv(res, file = result_name) 
  results_list[[result_name]] <- res 
} 
