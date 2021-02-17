# New file doing WGCNA analysis
#
BiocManager::install("WGCNA")

library(stringr)
library(WGCNA)
library(dplyr)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()


# load data ---------------------------------------------------------------
# Males
#norm_cts <- read.csv("input_data/normalized_counts_dseq_males.csv",
#                     row.names = 1) %>% as.data.frame %>% t

norm_cts <- read.csv("input_data/normalized_counts_dseq_all.csv",
                     row.names = 1) %>% as.data.frame %>% t

# Extract phenotype data
rownames(norm_cts)
phenotypes <- data.frame(row.names = row.names(norm_cts), sample =row.names(norm_cts) ) %>%
        mutate(
                species = str_sub(sample,-3),
                treatment = str_sub(sample, start = -4, end = -4),
                sp_tret = str_sub(sample,-4))

# Subset data just to fortis
fortis <- phenotypes$sample[phenotypes$species=="FOR"]
cra <- phenotypes$sample[phenotypes$species=="CRA"]
norm_cts <- norm_cts[rownames(norm_cts) %in% cra,]
dim(norm_cts)

# Inspect and compile data  -----------------------------------------------
# Check to see if there are genes that have 0 variance
# Returns a logical vector T/F for each gene, and T/F for each sample
gsg <- goodSamplesGenes(norm_cts, verbose = 3)
gsg$allOK # is everything ok?

# Remove bad genes; sample rows with genes passing QC and samples passing QC
norm_cts <- norm_cts[gsg$goodSamples, gsg$goodGenes]


# Gene expression network construction ------------------------------------

powers = c(c(1:10), seq(from = 12, to=20, by=2))


sft = pickSoftThreshold(norm_cts, powerVector = powers, verbose = 5)
# Plot the results:
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Looks like 12? is our best soft threshold power

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

# Change this accordingly
softPower = 12
adjacency <- adjacency(norm_cts, power = softPower)

# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")


# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

#This sets the minimum number of genes to cluster into a module
minModuleSize = 30

# Dynamic tree cutting, which clusters tree into modules
dynamicMods = cutreeDynamic(
        dendro = geneTree,
        distM = dissTOM,
        deepSplit = 2,
        pamRespectsDendro = FALSE,
        minClusterSize = minModuleSize
)
table(dynamicMods) # look at number of genes per module


dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(norm_cts, colors = dynamicColors, softPower = softPower)
MEs = MEList$eigengenes
MEDiss = 1 - cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
#save(dynamicMods, MEList, MEs, MEDiss, METree, file = "Network_allSamples_signed_RLDfiltered.RData")

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(norm_cts, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

pdf(file = "output_plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors
# Associate clusters with phenotypes --------------------------------------

nGenes = ncol(norm_cts)
nSamples = nrow(norm_cts)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(norm_cts, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Make phenotypes numeric
num_pheno <- phenotypes %>% filter(species == "FOR") %>% dplyr::select(treatment) %>%
        mutate_all(as.factor) %>% mutate_all(as.numeric)


moduleTraitCor = cor(MEs, num_pheno, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
colnames(norm_cts)[moduleColors == "darkorange2"] %>% write.csv("output_plots/wgcna_cra_darkorange_modgenes.csv")


