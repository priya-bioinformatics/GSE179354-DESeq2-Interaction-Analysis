# Set working directory
setwd("D:/Bioinformatics project/DESeq2/Counts")

getwd()
list.files()
files <- list.files(pattern = "counts.txt")
length(files)
#Reads one file to verify its structure
test <- read.table(files[1], header = FALSE)
head(test)
#Function to read one count file
read_counts <- function(file) {
  df <- read.table(file, header = FALSE)
  colnames(df) <- c("GeneID", gsub(".counts.txt", "", file))
  df
}
#Read all samples at once
count_list <- lapply(files, read_counts)
#Merge all samples into one matrix
count_matrix <- Reduce(
  function(x, y) merge(x, y, by = "GeneID"),
  count_list
)
#Prepare matrix for DESeq2
rownames(count_matrix) <- count_matrix$GeneID
count_matrix <- count_matrix[, -1]
#Check dimensions
dim(count_matrix)
#Filter low-expression genes
keep <- rowSums(count_matrix >= 10) >= 4
count_matrix_filt <- count_matrix[keep, ]
#Load Metadata
install.packages("readxl")
library(readxl)
colData <- read_excel("Metadata.xlsx")
#Prepare Metadata for DESeq2
colData <- as.data.frame(colData)
rownames(colData) <- colData$sample
colData$sample <- NULL
#Confirm Sample order matches
all(colnames(count_matrix_filt) == rownames(colData))
#Convert variables to factors
colData$cytokine <- factor(colData$cytokine)
colData$FP <- factor(colData$FP)
str(colData)
#Load DESeq2
library(DESeq2)
#Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filt,
  colData = colData,
  design = ~ FP + cytokine
)
#Run DESeq2 pipeline
dds <- DESeq(dds)
# TNFa vs Control
res_TNFa <- results(dds, contrast = c("cytokine", "TNFa", "Control"))
summary(res_TNFa)
#IFNg vs Control
res_IFNg <- results(dds, contrast = c("cytokine", "IFNg", "Control"))
summary(res_IFNg)
#TNFa + IFNg vs Control
res_combo <- results(dds, contrast = c("cytokine", "TNFa_IFNg", "Control"))
summary(res_combo)
#Save DE results
write.csv(as.data.frame(res_TNFa), "DE_TNFa_vs_Control.csv")
write.csv(as.data.frame(res_IFNg), "DE_IFNg_vs_Control.csv")
write.csv(as.data.frame(res_combo), "DE_TNFa_IFNg_vs_Control.csv")
# Create the DESeq2 dataset with interaction model
dds_int <- DESeqDataSetFromMatrix(
  countData = count_matrix_filt,
  colData = colData,
  design = ~ cytokine * FP 
)
# Run DESeq2
dds_int <- DESeq(dds_int)
# Interaction results for TNFa
res_FP_TNFa <- results(dds_int, name = "cytokineTNFa.FPNoFP")
summary(res_FP_TNFa)
# Check all coefficients
resultsNames(dds_int)
#Interaction effect for IFNg
res_FP_IFNg <- results(dds_int, name = "cytokineIFNg.FPNoFP")
summary(res_FP_IFNg)
#Interaction effect for TNFa+IFNg combination
res_FP_combo <- results(dds_int, name = "cytokineTNFa_IFNg.FPNoFP")
summary(res_FP_combo)
#Main FP effect at control level
res_FP_control <- results(dds_int, contrast = c("FP", "FP", "NoFP"))
summary(res_FP_control)
#Save interaction results as CSV
write.csv(as.data.frame(res_FP_TNFa), "FP_vs_TNFa.csv")
write.csv(as.data.frame(res_FP_IFNg), "FP_vs_IFNg.csv")
write.csv(as.data.frame(res_FP_combo), "FP_vs_TNFa_IFNg.csv")
write.csv(as.data.frame(res_FP_control), "FP_vs_Control.csv")
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
# List of results objects and contrast names
res_list <- list(
  "TNFa vs Control" = res_TNFa,
  "IFNg vs Control" = res_IFNg,
  "TNFa+IFNg vs Control" = res_combo,
  "FP vs Control" = res_FP_control,
  "FP effect on TNFa" = res_FP_TNFa,
  "FP effect on IFNg" = res_FP_IFNg,
  "FP effect on TNFa+IFNg" = res_FP_combo
)
# Function to count Up and Down genes
count_DE_genes <- function(res) {
  up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  return(c(Up = up, Down = down))
}

# Apply function to all results
summary_df <- lapply(res_list, count_DE_genes) %>% 
  bind_rows(.id = "Contrast")

# Reshape for ggplot
summary_long <- summary_df %>%
  pivot_longer(cols = c("Up", "Down"),
               names_to = "Regulation",
               values_to = "GeneCount")

# Plot
ggplot(summary_long, aes(x = Contrast, y = GeneCount, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Differentially Expressed Genes",
       x = "Contrast",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up" = "#E69F00", "Down" = "#56B4E9")) +
  theme_minimal(base_size = 14)
#Save the plot
ggsave(
  filename = " D:/Bioinformatics project/DESeq2/Counts /DEG_summary_barplot.pdf",
  width = 10,
  height = 6,
  dpi = 300
)
