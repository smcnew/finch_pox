# ok try cluster profiler

BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
install.packages("wordcloud")
install.packages("ggupset")

library(clusterProfiler)
library(wordcloud)
library(ggplot2)
library(DOSE)
library(enrichplot)
library(dplyr)
library(ggupset)

# data --------------------------------------------------------------------


# SET THE DESIRED ORGANISM HERE (chicken)
organism = "org.Gg.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

genes <- read.csv("input_data/for_i_u_clusterprofiler.csv")
original_gene_list <- genes$log2.FC.
names(original_gene_list) <- genes$GeneID
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

keytypes(org.Gg.eg.db)

# Gene set enrichment

eg <- bitr(genes$GeneID, fromType = "SYMBOL", toType = "EVIDENCE", OrgDb = organism)

gse <- gseGO(geneList = gene_list,
             ont = "BP",
             keyType = "SYMBOL",
             nPerm = 100,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "none"
             )

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

# GO analysis
# Pull out "significant" genes and see what they're enriched for
# an example
# enrichGO(sample_gene, OrgDb=org.Gg.eg.db, pvalueCutoff=1, qvalueCutoff=1) %>% head
sig_dif <- genes$GeneID[genes$significant=="yes"]

ego <- enrichGO(
  sig_dif,
  ont = "BP",
  OrgDb = org.Gg.eg.db,
  #readable = TRUE,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  keyType = "SYMBOL"
)
as.data.frame(ego) %>% View

ego <- setReadable(ego, OrgDb = org.Gg.eg.db)



# various visualizations
upsetplot(ego)

barplot(ego,
        drop = TRUE,
        showCategory = 10,
        title = "GO Biological Pathways",
        font.size = 8)
dotplot(ego)
emapplot(ego)
goplot(ego, showCategory = 10, max.overlaps = 3)

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(ego, categorySize="pvalue", foldChange=gene_list)

# KEGG Enrichment
sig_dif_trans <- bitr(sig_dif, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = organism)
dim(sig_dif_trans)
keytypes(org.Gg.eg.db)

# Translate whole list to ENTREZID and remove duplicates
# Starts with original_gene_list: a named vector of log2fold changes

ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "UNIPROT", OrgDb=organism) # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
head(dedup_ids)
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = genes[genes$GeneID %in% dedup_ids$SYMBOL,]
colnames(df2)[1] <- "SYMBOL"
df2 <- merge(df2, dedup_ids) #add kegg gene name

kegg_gene_list <- df2$log2.FC.  #list of known genes
names(kegg_gene_list) <- df2$UNIPROT
kegg_gene_list <- na.omit(kegg_gene_list) %>% sort(decreasing = https://david.ncifcrf.gov/helps/2D_Introduction_files/conversion_result.jpgTRUE) # get rid of NAs

kegg_genes <- filter(df2, significant == "yes") %>% pull(UNIPROT)


# kegg code is gfr
kegg_organism = "gfr"
kk <- enrichKEGG(
  gene = kegg_genes,
  universe = names(kegg_gene_list),
  organism = "gfr",https://david.ncifcrf.gov/helps/2D_Introduction_files/conversion_result.jpg
  pvalueCutoff = 0.05,
  keyType = "uniprot" #one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
)
head(genes)

