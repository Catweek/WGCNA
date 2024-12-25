# Set the working directory to the location of the input files
setwd("Yourfolder/WGCNA_QC")

#install.packages("WGCNA")            # Install in case you don't have these and load necessary libraries
library(WGCNA)
num_threads <- 4                      # Set the number of threads to, for example, 4
allowWGCNAThreads(nThreads = num_threads)
options(stringsAsFactors = FALSE)     # Set options to prevent strings from being converted to factors

# NOTE: Uncomment this block if your input data needs normalization; # Uncomment and edit the following lines if needed
# install.packages("DESeq2")
# library(DESeq2) # Required for DESeqDataSetFromMatrix and normalization
# dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data0[,-181], colData = sample_metadata, design = ~ Zone)
# mcols(dataExpr_deseq)$basepairs = data0$geneLengt1
# fpkm_matrix <- fpkm(dataExpr_deseq)

# Step 1: Read the gene counts table
cat("Reading gene counts table...\n")
data0 <- read.table("gene_counts_table_WGCNA.txt", header = TRUE, row.names = 1, sep = "\t")

# Display the structure and dimensions of the gene counts table
cat("Gene counts table dimensions: ", dim(data0), "\n")
cat("First 5 rows of the gene counts table:\n")
print(head(data0, 5))
datExpr <- t(data0)    # Transpose the data so that samples are in rows and genes are in columns
print(datExpr[1:5, 1:5])

# Load sample metadata (ensure sample_info.csv contains the `sample_ID` column)
cat("Reading sample metadata...\n")
sample_metadata <- read.csv(file = "sample_info.csv")
print(head(sample_metadata, 5))   # Display the structure and a preview of the sample metadata

# Check if sample IDs in metadata match column names in the gene counts table
cat("Matching sample IDs in metadata with column names in the gene counts table...\n")
id_match <- match(sample_metadata$sample_ID, colnames(data0))
if (any(is.na(id_match))) {
  cat("Warning: Some sample IDs in metadata do not match column names in the gene counts table.\n")
} else {
  cat("All sample IDs in metadata match column names in the gene counts table.\n")
}
# datExpr <- datExpr[, 1:5000]  # NOTE: Uncomment this line to reduce the dataset to a manageable number of genes for testing
# Step 2: Dendrogram: to calculate sample distances and cluster the samples
cat("Calculating sample distances and performing hierarchical clustering...\n")
sampleTree <- hclust(dist(datExpr), method = "average")
sizeGrWindow(12, 9)               #Step 3: Visualize the sample clustering to detect outliers
par(cex = 0.6)  # Set the size of plotted text
par(mar = c(0, 4, 2, 0))  # Set the margin sizes
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

# Choose a set of soft threshold parameters | Scale-free topology fit index as a function of the soft-threshold power
powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

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
# unsigned -> nodes with positive & negative correlation are treated equally 
# signed -> nodes with negative correlation are considered *unconnected*, treated as zero

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
pdf(file = "4-module_tree_blockwise_30.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

################################################################################
################################################################################
# Option 2a: step-by-step
power = power
adjacency = adjacency(datExpr, power = power)
TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
dissTOM = 1-TOM


# Option 2b: 
TOM = TOMsimilarityFromExpr(datExpr, power = power)
dissTOM = 1-TOM 
dim(dissTOM)

#===============================================================================
#
#  Construct modules 
#
#===============================================================================
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# Module identification using dynamic tree cut
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30);
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#  Merge modules
#
#===============================================================================
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Merge close modules
MEDissThres=0.40
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
pdf(file = "5-merged_Module_Tree_again.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()
write.table(merge$oldMEs,file="oldMEs1.txt");
write.table(merge$newMEs,file="newMEs1.txt");

# Read clinical data as traits
bac_traits = read.table("clinical_traits.txt", header = T, sep = "\t")
rownames(bac_traits) = bac_traits[, 1]
bac_traits = bac_traits[, -1]
# sample names should be consistent in eigen genes and traits !!!!
bac_traits = bac_traits[match(rownames(MEs), rownames(bac_traits)), ]
table(rownames(MEs) == rownames(bac_traits))

MEs_path <- "severe.txt"
MEs_new <- read.table(MEs_path, header = TRUE)
# Make sure sample names are consistent between MEs and traits
bac_traits <- bac_traits[match(rownames(MEs_new), rownames(bac_traits)), ]
table(rownames(MEs_new) == rownames(bac_traits))

# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs_new, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.table(moduleTraitCor,file="moduleTrait_correlation_11.txt");
write.table(moduleTraitPvalue,file="moduleTrait_pValue_11.txt");


#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("sev11.pdf", width = 35, height = 30)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
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

############# Summary ###################################

head(datExpr)[1:5,1:5] # transcriptome data

head(sample_metadata)[1:5,] # metadata (sample info)
head(bac_traits)[1:5,1:5] # external trait


#=====================================================================================
#
#   Cytoscape
#
#=====================================================================================


if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}

#https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/
library(RCy3)

cytoscapePing () # make sure cytoscape is open
cytoscapeVersionInfo ()

###### for yellow module of the merged data (newMEs) #################################
edge <- read.delim("output_for_cytoscape_new/merge_CytoscapeInput-edges-grey.txt")
colnames(edge)
colnames(edge) <- c("source", "target","weight","direction","fromAltName","toAltName")

node <- read.delim("output_for_cytoscape_new/merge_CytoscapeInput-nodes-grey.txt")
colnames(node)  
colnames(node) <- c("id","altName","node_attributes") 

createNetworkFromDataFrames(node,edge[1:50,], title="my first network", collection="DataFrame Example")

################ customise the network visualization ##################################
# use other pre-set visual style
setVisualStyle('Marquee')

# set up my own style
style.name = "myStyle"
defaults <- list(NODE_SHAPE="diamond",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','node_attributes','d',c("A","B"), c("#FF9900","#66AAAA"))
arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
edgeWidth <- mapVisualProperty('edge width','weight','p')

createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes,edgeWidth))
setVisualStyle(style.name)