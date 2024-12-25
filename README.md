# WGCNA
This repository contains R scripts and associated data used for performing Weighted Gene Co-expression Network Analysis (WGCNA). The analysis is aimed at exploring gene expression patterns across samples, identifying modules of highly correlated genes, and examining their relationships with clinical traits or other factors.

This repository has input files and R script needed. 
All that is required for the tutorial in the Youtube Channel (Liquid Brain, https://www.youtube.com/c/LiquidBrain)

It covers: 
1. What data you need for WGCNA
2. How to perform network construction and module detection
3. Correlate the modules with external trait (continuous type)
4. Further investigation on particular module-trait relationship
5. Visualization (e.g. heatmap)
6. Summary of the whole WGCNA analysis

## Data Overview
### Data Normalization (Must)
My data was normalised. However if your data is not normalsied, you must do it with DESeq2 (codes are attached in .R). Just uncomment and modify to use the `DESeq2` package. Normalization is important to account for differences in library sizes or other technical biases in RNA-seq data. 
### Gene Counts Table
The primary input data used in this analysis is a gene expression count table, which contains gene expression data from RNA-seq experiments. 
The .txt table is structured as follows: 
![image](https://github.com/user-attachments/assets/e8d1bec6-ffcf-4656-9fa6-9db00882ee59)

### Step 1: Reading the Gene Counts Table
After loading the data, we transpose it using `t(data0)` so that the samples are arranged as rows and genes as columns. This transposition is required for subsequent WGCNA analysis, where genes are represented as columns and samples as rows.
```r
cat("Reading gene counts table...\n")
data0 <- read.table("gene_counts_table_WGCNA.txt", header = TRUE, row.names = 1, sep = "\t")
cat("Gene counts table dimensions: ", dim(data0), "\n")
cat("First 5 rows of the gene counts table:\n")
print(head(data0, 5))
datExpr <- t(data0)    # Transpose the data so that samples are in rows and genes are in columns
print(datExpr[1:5, 1:5])
```
The first 5 rows of the gene counts table and a small sample of the transposed data are printed for verification.

The analysis also uses sample metadata, which contains information about each sample (e.g., experimental conditions, clinical traits). In my case i had disease severity (referred as Zones: Z1, Z2, Z3 as mild, moderate and severe) 
```r
cat("Reading sample metadata...\n")
sample_metadata <- read.csv(file = "sample_info.csv")
print(head(sample_metadata, 5))   # Display the structure and a preview of the sample metadata
```
![image](https://github.com/user-attachments/assets/46b62ec4-9427-4d57-8b35-d5bb406396a8)

#### Now who wants chaos in analysis? 
Not me! So lets just verifyy things between both of our input files. The script performs a check using the `match()` function to align the `sample_ID` column from the metadata with the column names in the gene counts table (the one before transposing - data0). 
```r
cat("Matching sample IDs in metadata with column names in the gene counts table...\n")
id_match <- match(sample_metadata$sample_ID, colnames(data0))
if (any(is.na(id_match))) {
  cat("Warning: Some sample IDs in metadata do not match column names in the gene counts table.\n")
} else {
  cat("All sample IDs in metadata match column names in the gene counts table.\n")
}
```
If there is a mismatch, it indicates an issue with the file formatting/sample labeling, and will need to be addressed before proceeding with analysis.
To speed up the analysis or for testing purposes, the dataset can be reduced to a smaller set of genes by uncommenting the following line. This will keep only the first 5000 genes for analysis: I will let it be uncommented onyl because obviously i am having human data with like 20K genes approx.
```r
# datExpr <- datExpr[, 1:5000]  # NOTE: Uncomment this line to reduce the dataset to a manageable number of genes for testing
```
### Step 2: Dendrogram – Calculating Sample Distances and Clustering
Now, I calculate the distances between samples based on their gene expression profiles and perform hierarchical clustering to identify potential outliers. This is crucial for assessing the overall similarity between samples and detecting any samples that deviate significantly from the others.

The `hclust()` function is used to perform hierarchical clustering on the distance matrix (`dist(datExpr)`). 
```r
cat("Calculating sample distances and performing hierarchical clustering...\n")
sampleTree <- hclust(dist(datExpr), method = "average")
```
The `plot()` function is used to generate the plot, with red horizontal lines indicating the thresholds for outlier detection.
```r
sizeGrWindow(12, 9)               # Set the size of the plot window
par(cex = 0.6)                    # Set the size of plotted text
par(mar = c(0, 4, 2, 0))          # Set margin sizes
plot(sampleTree, 
     main = "Dendrogram to look for outliers", 
     sub = "", 
     xlab = "", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     cex.main = 2)
#abline(h = 130, col = "red")
#abline(h = 150, col = "red")
#abline(h = 190, col = "red")
```
![w3](https://github.com/user-attachments/assets/95a1f513-b8d6-4b24-9b5d-028163bc4fbe)

The optional red lines in the dendrogram indicate where samples can be split to detect outliers. Adjust the cutoff as needed depending on your dataset.

### Step 3: Soft Thresholding for Network Construction
In WGCNA, choosing the right soft threshold is critical for constructing a network that captures meaningful gene co-expression relationships while maintaining a scale-free topology. This is where the above three red lines help you. The `pickSoftThreshold()` function is used to determine the optimal soft threshold based on the scale-free topology fit index.

#### Picking the Soft Threshold
A range of soft threshold powers is tested to select the optimal value for network construction - and i did it like times minimum until i found the ultimate automatic option. The soft threshold is chosen by analyzing the scale-free topology fit index (R²) and ensuring that it meets a minimum threshold, typically 0.90:
```r
powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
```
Two plots are generated: one for the scale-free topology fit index and another for the mean connectivity across different soft thresholds.
```r
pdf(file = "2-n-sft.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
```
![image](https://github.com/user-attachments/assets/19731ac3-5ae3-4c54-9467-9270521149c0)
These plots help determine the optimal power for the soft threshold, ensuring that the resulting network is both scale-free and well-connected.

### Step 4: Constructing the Network
Once the soft threshold is selected, we use the `blockwiseModules()` function to construct the network and identify gene modules. The function takes the expression data, the selected power, and various other parameters for network construction:
```r
# Turn data expression into topological overlap matrix
power=sft$powerEstimate #4

# Option 1: automatic
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = power,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)
cor<- stats::cor
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
pdf(file = "4-module_tree_blockwise_30.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
```
The `blockwiseModules()` function divides genes into modules based on their co-expression patterns.

### Step 5: Topological Overlap Matrix (TOM)
For further analysis, the topological overlap matrix (TOM) is computed. TOM measures the similarity between genes based on their co-expression patterns, which is essential for identifying modules of co-expressed genes.

#### Option 2a: Automatic TOM Calculation
In this approach, TOM is calculated directly from the adjacency matrix using the `TOMsimilarity()` function:
```r
power = power
adjacency = adjacency(datExpr, power = power)
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM
```
#### Option 2b: Manual TOM Calculation
Alternatively, TOM can be computed directly from the expression data using the `TOMsimilarityFromExpr()` function:
```r
TOM = TOMsimilarityFromExpr(datExpr, power = power)
dissTOM = 1 - TOM
dim(dissTOM)
```
The resulting `dissTOM` matrix represents the dissimilarity between genes, which can be used for further network analysis.

## Step 6: Construct Modules

### Module Identification Using Dynamic Tree Cut
First, we calculate the gene tree using hierarchical clustering on the topological overlap matrix (TOM), and then we plot the gene dendrogram.
We then identify gene modules using a dynamic tree cut method. In this case, we set the minimum module size to 30 to ensure that only large modules are considered.
```r
# Module identification using dynamic tree cut
geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, 
                            pamRespectsDendro = FALSE, minClusterSize = 30)
table(dynamicMods)
length(table(dynamicMods))
```
Next, we convert the numeric module labels into colors and plot the dendrogram with module colors.
```r
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file = "4-module_tree.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()
```
![image](https://github.com/user-attachments/assets/0afa1016-6d70-4753-aab1-5a5dc838693e)

## Step 7: Merge Modules

### Calculate Module Eigengenes
Next, we calculate the eigengenes of the modules, which summarize the gene expression patterns within each module.
```r
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
```
### Module Dissimilarity and Clustering
We compute the dissimilarity between module eigengenes and cluster them to assess how similar the modules are to each other.
```r
# Calculate dissimilarity of module eigengenes
MEDiss = 1 - cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
```
### Merge Close Modules
We merge modules that are highly correlated (i.e., those with a module eigengene dissimilarity less than a threshold of 0.40).
```r
# Merge close modules
MEDissThres = 0.40
abline(h = MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
# Plot merged module tree
pdf(file = "5-merged_Module_Tree_again.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Save the old and new module eigengenes
write.table(merge$oldMEs, file = "oldMEs1.txt")  before merging scores
write.table(merge$newMEs, file = "newMEs1.txt")  after merging scores
```
![image](https://github.com/user-attachments/assets/d8547c22-1975-4698-a1bf-f04f0b4ed601)

### Calculate Pearson Correlations Between Module Eigen-genes and Clinical Traits
![image](https://github.com/user-attachments/assets/4cf6428e-8361-48a7-ad8a-925ae48a7e91)

bac_traits = read.table("clinical_traits.txt", header = T, sep = "\t")
rownames(bac_traits) = bac_traits[, 1]
bac_traits = bac_traits[, -1]
### sample names should be consistent in eigen genes and traits !!!!
bac_traits = bac_traits[match(rownames(MEs), rownames(bac_traits)), ]
table(rownames(MEs) == rownames(bac_traits))
MEs_path <- "severe.txt"     #for one severity per time as we are correlating each disease severity with clinical parameter
MEs_new <- read.table(MEs_path, header = TRUE)
![image](https://github.com/user-attachments/assets/2088bbc9-924b-4ade-9a1e-6feac02d7009)

### Make sure sample names are consistent between MEs and traits
bac_traits <- bac_traits[match(rownames(MEs_new), rownames(bac_traits)), ]
table(rownames(MEs_new) == rownames(bac_traits))
### PCC
moduleTraitCor = cor(MEs_new, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.table(moduleTraitCor,file="moduleTrait_correlation_11.txt");
write.table(moduleTraitPvalue,file="moduleTrait_pValue_11.txt");
### Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
### Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("plot.pdf", width = 35, height = 30)
par(mar = c(15, 12, 5, 5));
### Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(bac_traits),
               yLabels = colnames(MEs_new),
               ySymbols = colnames(MEs_new),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.8,
               cex.axis = 5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

![image](https://github.com/user-attachments/assets/896b0955-1735-448b-b5e2-0d6896b1f218)

Severe data association with clincial parameters

### Summary:
head(datExpr)[1:5,1:5] # transcriptome data
head(sample_metadata)[1:5,] # metadata (sample info)
head(bac_traits)[1:5,1:5] # external trait

### References:

The original tutorial provided by the creators (Peter Langfelder and Steve Horvath) https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
A nice Nature Plants paper by Yu et al. https://www.nature.com/articles/s41477-021-00897-y (Check their "data availability" section for the link to the github R script)
