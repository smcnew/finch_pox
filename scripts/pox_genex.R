# A script to estimate expression counts from alignment data

#BiocManager::install("DESeq2")
#BiocManager::install("biomaRt")
#BiocManager::install('EnhancedVolcano')

library(DESeq2)
library(biomaRt)
library(stringr)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
# Data --------------------------------------------------------------------

# TREX produced file of raw counts
raw_counts <- read.table("TREX_Results/2008D-Gfortis_rawCounts.txt")
head(raw_counts)

# idk what's happening here
#sig_diff <- read.csv("input_data/sig_diff.csv")
#head(sig_diff)
#genes <- sig_diff %>% pull %>% as.character()

# for vs. cra genes
ifor_cra <- read.csv("TREX_Results/ifor_vs_icra.csv")
ufor_cra <- read.csv("TREX_Results/ufor_vs_ucra.csv")
# Convert gene to ensembl -------------------------------------------------
# what is this??

ensembl = useMart("ensembl", dataset = "tguttata_gene_ensembl")

getBM(attributes =
      filters = 'ensembl_gene_id',
      values = genes,
      mart = ensembl)
listAttributes(ensembl) %>% dplyr::filter(., name =="hgnc_symbol")


# DESeq2 ------------------------------------------------------------------
# All samples (males and females)
# Normalize counts for analysis in wgcna

raw_counts

# Pull metadata out
metadata <- data.frame(sample = colnames(raw_counts)) %>%
  dplyr::mutate(
    species = str_sub(sample,-3),
    treatment = str_sub(sample, start = -4, end = -4),
    sp_tret = str_sub(sample,-4),
    band = unlist(str_split(sample, "[.]"))[2]
  )


# need to make sure metadata and columns are in same order
all(metadata$sample == colnames(raw_counts)) #TRUE
dds_all <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = metadata,
                              design = ~ sp_tret)

# Remove genes with fewer than 10 reads
keep <- rowSums(counts(dds_all)) >= 10
dds_all <- dds_all[keep,]
dds_all <- DESeq(dds_all)

#
res_for <- results(dds_all, contrast = c("sp_tret", "IFOR", "UFOR"), alpha = .05)
summary(res_for)

res_cra <- results(dds_all, contrast = c("sp_tret", "ICRA", "UCRA"), alpha = .05)
summary(res_cra)

res_cra_for <- results(dds_all, contrast = c("sp_tret", "IFOR", "ICRA"), alpha = .05)
summary(res_cra_for)

# Normalize reads and export
vst(dds_all, blind = TRUE) %>% assay %>% write.csv("input_data/normalized_counts_dseq_all.csv")
vsd_all <- vst(dds_all, blind = TRUE) #create this object for PCA

# Volcano plots: all individuals  --------------------------------------------------


pdf("output_plots/volcanos_m_andf.pdf")
par(mfrow=c(1,2))
EnhancedVolcano(res_for,
                lab = rownames(res_for),
                title = "Medium ground finch",
                x = 'log2FoldChange',
                y = 'pvalue',
                #pCutoff = .01,
                FCcutoff = 2,
                xlim = c(-7,7),
                col=c('black', 'black', 'blue', 'red3'),#'#009988'
                colAlpha = 0.5,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = .8)


EnhancedVolcano(res_cra,
                lab = rownames(res_cra),
                title = "Vegetarian finch",
                x = 'log2FoldChange',
                y = 'pvalue',
                #pCutoff = .01,
                FCcutoff = 2,
                xlim = c(-7,7),
                col=c('black', 'black', 'blue', 'red3'),#'#009988'
                colAlpha = 0.5,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = .8)

dev.off()




# Males only --------------------------------------------------------------
males_raw <- dplyr::select(raw_counts, contains(".M."))
males_metadata <- data.frame(sample = colnames(males_raw)) %>%
  mutate(species = str_sub(sample, -3),
         treatment = str_sub(sample, start = -4, end = -4),
         sp_tret = str_sub(sample, -4),
         band = unlist(str_split(sample, "[.]"))[2])

# need to make sure metadata and columns are in same order
all(males_metadata$sample == colnames(males_raw)) #TRUE
head(males_metadata)
dds <- DESeqDataSetFromMatrix(countData = males_raw,
                              colData = males_metadata,
                              design = ~ sp_tret)
# Pre-filter to remove low count genes
keep <- rowSums(counts(dds)) >= 10 #clean data by removing rows with < 10 cts
dds <- dds[keep,]
dds <- DESeq(dds)
res_for <- results(dds, contrast = c("sp_tret", "UFOR", "IFOR"), alpha = .05)
res_cra <- results(dds, contrast = c("sp_tret","UCRA", "ICRA"), alpha = .05)

# Then use assay to extract the matrix of normalized counts
vsd <- vst(dds, blind = T)
assay(vsd) %>% write.csv("input_data/normalized_counts_dseq_males.csv")


# PCA -------------------------------------------------


pca <- plotPCA(vsd_all, intgroup = c("species", "treatment"), returnData = T)
percentVar <- round(100 * attr(pca, "percentVar"))
ggplot(pca, aes(PC1, PC2, color=species, shape = treatment)) +
  geom_point(size = 3, stroke = 2, alpha = 0.9) +
  scale_color_manual(values = c("#33BBEE", "#EE3377"), labels = c("Vegetarian","Med. Ground")) +
  scale_shape_manual(values = c(3, 1), labels = c("Infected", "Uninfected")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance"))+
  ylab(paste0("PC2: ", percentVar[2], "% variance"))+
  theme(axis.text=element_text(size=18),
  axis.title=element_text(size=20))

ggsave("output_plots/pca.pdf", width = 7, height = 5)




#

