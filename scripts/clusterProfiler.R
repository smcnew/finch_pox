# ok try cluster profiler

# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
# install.packages("wordcloud")
# install.packages("ggupset")
# BiocManager::install("mygene")

library(clusterProfiler)
library(wordcloud)
library(ggplot2)
library(ggridges)
library(DOSE)
library(enrichplot)
library(ggupset)
library(magrittr)
library(tibble)
library(mygene) #to translate gene names among types
library(ggnewscale)
library(msigdbr) #for msigdb databases
library(stringr)
library(cowplot)
library(VennDiagram)
library(dplyr)

set.seed(1987)
# data --------------------------------------------------------------------


# Results from DESeq2, filtered from TREX results
# Sample sizes:
# CRA: 6 infected M, 9 uninfected M;
# FOR: 6 infected M, 7 uninfected M; 4 infected F, 3 uninfected F
#
# Males and females separated
for_miu <- read.csv("input_data/for_m_i_u_clusterprofiler.csv") # fortis infected vs not
cra_miu <- read.csv("input_data/cra_m_i_u_clusterprofiler.csv") # cra i vs u
for_fiu <- read.csv("input_data/for_f_i_u_clusterprofiler.csv") # female for i vs u
for_cra_mi <- read.csv("input_data/for_i_cra_i_m_clusterprofiler.csv") %>%
  mutate(log2.FC. = log2.FC.*(-1)) # Cra I vs For I; multiply by -1 to make Cra the reference
for_cra_mu <- read.csv("input_data/for_u_cra_u_m_clusterprofiler.csv") %>%
  mutate(log2.FC. = log2.FC.*(-1)) # Cra I vs For I; cra is ref.

# list all the items
list_sex <- c("for_miu", "cra_miu", "for_fiu", "for_cra_mi", "for_cra_mu")

# males and females combined

for_iu <- read.csv("input_data/for_i_u_clusterprofiler.csv")
cra_iu <- read.csv("input_data/cra_i_u_clusterprofiler.csv")
cra_for_i <- read.csv("input_data/cra_i_for_i_clusterprofiler.csv") %>%
mutate(log2.FC. = log2.FC.*(-1)) # Cra I vs For I; cra is ref; multiply by -1 to make Cra the ref.
# i.e. "Compared to veggie, gound were up-regulated"

cra_for_u <- read.csv("input_data/cra_u_for_u_clusterprofiler.csv") %>%
  mutate(log2.FC. = log2.FC.*(-1)) # Cra I vs For I; cra is ref; multiply by -1 to change to For is ref

# list all the items
list_all <- c("for_iu", "cra_iu", "cra_for_i", "cra_for_u")


# summary stats -----------------------------------------------------------
# how many genes were significant?
for(i in list_sex) {
 # print(i)
get(i) %>%
    pull(DIFF_EXP) %>% sum %>% paste(i, .) %>% print
}

for(i in list_all) {
  # print(i)
  get(i)  %>% filter(sex.linked =="-") %>%
pull(DIFF_EXP) %>% sum %>% paste(i, .) %>% print
}



# SET THE DESIRED REF ORGANISM HERE (chicken); not sure I need this anymore.
#organism = "org.Gg.eg.db"
#BiocManager::install(organism, character.only = TRUE)
#library(organism, character.only = TRUE)

# if library doesn't load try
# options(connectionObserver = NULL) then
# call dbDisconnect() when finished working with a connection (this didn't work)

# functions ---------------------------------------------------------------
# Need functions to translate gene names and return either a) sorted gene list
# or b) decreasing sorted vector (<- for GSEA)


# Finch IDs
ids <- queryMany(as.character(for_iu$GeneID), scope = "symbol", fields = "entrezgene", species = "48883") %>%
  as.data.frame %>% dplyr::select(query, entrezgene) %>% filter(!is.na(entrezgene)) %>%
  dplyr::rename(GeneID = query)

# Chicken IDs
ids_gallus <- queryMany(as.character(for_iu$GeneID), scope = "symbol", fields = "entrezgene", species = "9031") %>%
  as.data.frame %>% dplyr::select(query, entrezgene) %>% dplyr::filter(!is.na(entrezgene)) %>%
  dplyr::rename(GeneID = query)

# sort_genes, returns a named vector of sorted fold changes, args = data and species
sort_genes <- function(data, species = "gallus") {
  if(species == "gallus") {
    idlist <- ids_gallus
  }
  else
    idlist <- ids
  dplyr::filter(data, GeneID %in% idlist$GeneID) %>%
    merge(idlist) %>%
    dplyr::select(log2.FC., entrezgene) %>%
    arrange(desc(log2.FC.)) %>% pull(log2.FC., name = entrezgene)
}
for_iu[for_iu$sex.linked=="-",] %>%
  dplyr::select(log2.FC., GeneID) %>% arrange(desc(log2.FC.)) %>% filter(GeneID == "HIF1A")

# pull_genes, returns a vector of entrez gene ids, args = data, species, log change.
pull_genes <- function (data, species = "gallus", log_change = 2) {
  if(species == "gallus") {
    idlist <- ids_gallus
  }
  else
    idlist <- ids
  filter(data, GeneID %in% idlist$GeneID) %>%
    merge(idlist) %>%
    filter(log2.FC. <= log_change) %>%
    pull(entrezgene)
}

# turn capitalized names into uncapitalized
uncap <- function(x) {
  str_replace(x, "HALLMARK_","") %>% str_replace_all ("_", " ") %>%
    str_to_lower
}

# MSigDB ------------------------------------------------------------------

# Download hallmark gene sets for chicken
h_gene_sets <- msigdbr(species = "Gallus gallus", category = "H")
#glimpse(h_gene_sets)

# Pull the hallmark set and the gene name
msigdbr_2gene <- h_gene_sets %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

# Download immune gene sets for chicken
h_immune_sets <- msigdbr(species = "Gallus gallus", category = "C7")

msigdbr_2gene_immune <- h_immune_sets %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

# Gene set enrichment with hallmark sets  ---------------------------------
#


# 1. fortis infected vs uninfected

for_iu_gsea <- sort_genes(for_iu[for_iu$sex.linked=="-",], species = "gallus") %>%
  GSEA(., TERM2GENE = msigdbr_2gene, pAdjustMethod = "fdr")
for_iu_gsea %>% dim # 14 sig dif
for_iu_gsea %>% ridgeplot(label_format = function(x) uncap(x))+
  ggtitle("Ground: Infected vs uninfected") +
  xlab(expression(log[2]*" fold change" ))

# fortis infected uninfected immune genes
#
#sort_genes(for_iu[for_iu$sex.linked=="-",], species = "gallus") %>%
#  GSEA(., TERM2GENE = msigdbr_2gene_immune, pAdjustMethod = "fdr") %>% ridgeplot()

# # looking into hypoxia a little more
# for_iu %>% filter(GeneID == "HIF1A")
# hypogenes <-
#   for_iu_gsea@result %>% filter(ID == "HALLMARK_HYPOXIA") %>%
#   pull(core_enrichment) %>% str_split(., "/") %>%
#   unlist %>% as.numeric
# hypogeneids <-
#   filter(ids_gallus, entrezgene %in% hypogenes) %>% pull(GeneID)
# for_iu %>% filter(GeneID %in% hypogeneids) %>%
#   arrange(desc(log2.FC.))
# dplyr::filter(data, GeneID %in% idlist$GeneID) %>%
#   merge(idlist) %>%
#   dplyr::select(log2.FC., entrezgene) %>%
#   arrange(desc(log2.FC.)) %>% pull(log2.FC., name = entrezgene)

#
#for_iu_gsea %>% gseaplot2(., geneSetID = 1:14)

# ## ggplot bar
# ggplot(for_iu_gsea@result, aes(reorder(ID, NES), NES)) +
#   geom_col(aes(fill=p.adjust<0.05)) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="Hallmark pathways NES from GSEA") +
#   theme_minimal()



# 2. fortis infected vs. cra infected
cra_for_i_gsea <- sort_genes(cra_for_i[cra_for_i$sex.linked=="-",], species = "gallus") %>%
  GSEA(., TERM2GENE = msigdbr_2gene, pAdjustMethod = "fdr")
ridgeplot(cra_for_i_gsea, label_format = function(x) uncap(x)) +  ggtitle('ground vs. veggie, all infected')
cra_for_i_gsea@result %>% dim # 9 sig diff
cra_for_i_gsea %>% gseaplot2(., geneSetID = 1:nrow(.))

# 3. cra infected vs uninfected
cra_iu_gsea <- sort_genes(cra_iu[cra_iu$sex.linked=="-",], species = "gallus") %>%
  GSEA(., TERM2GENE = msigdbr_2gene, pAdjustMethod = 'fdr')
cra_iu_gsea %>% ridgeplot(label_format = function(x) uncap(x))+  ggtitle('Veggie: infected vs. uninfected')
#cra_iu_gsea %>% gseaplot2(., geneSetID = 1:nrow(.))



# 4. fortis uninfected vs. cra uninfected
cra_for_u_gsea <- sort_genes(cra_for_u[cra_for_u$sex.linked=="-",], species = "gallus") %>%
  GSEA(., TERM2GENE = msigdbr_2gene, pAdjustMethod = 'fdr')
cra_for_u_gsea %>% ridgeplot(label_format = function(x) uncap(x)) +  ggtitle('uninfected fortis vs. uninfected veggie')
cra_for_u_gsea %>% gseaplot2(., geneSetID = 1:nrow(.))

# hallmark gsea just males  -----------------------------------------------

# Now males vs. females
# fortis infected vs. uninfected males
set.seed(1987)
for_miu_gsea <- sort_genes(for_miu, species = "gallus") %>%
  GSEA(., TERM2GENE = msigdbr_2gene, pAdjustMethod = "fdr")
for_miu_gsea
for_miu_gsea %>% ridgeplot() + ggtitle("ground infected vs. uninfected, males") #upsetplot
#for_miu_gsea %>% gseaplot2(., geneSetID = 1:nrow(.))

# fortis females infected vs uninfected
for_fiu_gsea <- sort_genes(for_fiu, species = "gallus") %>%
  GSEA(., TERM2GENE = msigdbr_2gene)
for_fiu_gsea@result[1:8,1:8]
for_fiu_gsea %>% ridgeplot() + ggtitle("fortis female i vs u")
for_fiu_gsea %>% gseaplot2(., geneSetID = 1:nrow(.))

# cra males infected vs uninfected.
cra_miu_gsea <- sort_genes(cra_miu, species = "gallus") %>%
  GSEA(., TERM2GENE = msigdbr_2gene)
ridgeplot(cra_miu_gsea) # nothing significant


# for males infected vs. cra males infected
for_cra_mi_gsea <- sort_genes(for_cra_mi) %>%
  GSEA(., TERM2GENE = msigdbr_2gene)
for_cra_mi_gsea %>% ridgeplot() + ggtitle("ground vs. veggie, infected males")

# fortis males uninfected vs. cra males uninfected
for_cra_mu_gsea <- sort_genes(for_cra_mu) %>%
  GSEA(., TERM2GENE = msigdbr_2gene)
for_cra_mu_gsea %>% ridgeplot()
for_cra_mu_gsea %>% gseaplot2(geneSetID = 1)



# hallmark ridgeplots for export ------------------------------------------

# all individuals
p1 <- for_iu_gsea %>% ridgeplot(label_format = function(x) uncap(x))+ ggtitle("Ground: infected vs uninfected") +
  xlab(expression(log[2]*" fold change" ))
p2 <- cra_iu_gsea %>% ridgeplot(label_format = function(x) uncap(x))+  ggtitle('Vegetarian: infected vs. uninfected')+
  xlab(expression(log[2]*" fold change" ))
p3 <- cra_for_i_gsea %>% ridgeplot(label_format = function(x) uncap(x)) +  ggtitle('Infected: ground vs. vegetarian')+
  xlab(expression(log[2]*" fold change" ))
p4 <- cra_for_u_gsea %>% ridgeplot(label_format = function(x) uncap(x)) +  ggtitle('Uninfected: ground vs. vegetarian')+
  xlab(expression(log[2]*" fold change" ))
pdf("output_plots/hallmark_ridgeplots.pdf", width = 12, height = 10)
plot_grid(p1, p2, p3, p4, labels = c('A', 'B', 'C', 'D'), label_size = 14)
dev.off()

# just males
# cra i vs. u no significant dif
p5 <- for_miu_gsea %>% ridgeplot(label_format = function(x) uncap(x))+ ggtitle("Ground: infected vs uninfected")+
  xlab(expression(log[2]*" fold change" ))
p6 <- for_cra_mi_gsea %>% ridgeplot(label_format = function(x) uncap(x)) +  ggtitle('Infected: ground vs. vegetarian')+
  xlab(expression(log[2]*" fold change" ))
p7 <- for_cra_mu_gsea %>% ridgeplot(label_format = function(x) uncap(x)) +  ggtitle('Uninfected: ground vs. vegetarian')+
  xlab(expression(log[2]*" fold change" ))
pdf("output_plots/hallmark_ridgeplots_males.pdf", width = 12, height = 10)
plot_grid(p5, p6, p7, labels = c('A', 'B', 'C'), label_size = 14)
dev.off()

# Venn diagram ------------------------------------------------------------

overlap <- calculate.overlap(x = list ("ground" = for_iu_gsea@result$ID,
                            "veggie" = cra_iu_gsea@result$ID,
                            "inf vs inv" = cra_for_i_gsea@result$ID))
draw.triple.venn(overlap)
vp <- venn.diagram(list (
                       "Uninf. vs. uninf.\nN = 3" = cra_for_u_gsea@result$ID,
                       "Inf. vs inf.\nN = 8" = cra_for_i_gsea@result$ID,
                       "Vegetarian\nN = 9" = cra_iu_gsea@result$ID,
                        "Ground\nN = 14" = for_iu_gsea@result$ID),
                   fill = 2:5, alpha = 0.3, cat.cex = 1.5, cex =1.4,
                   #print.mode = "percent", sigdigs = 2,
                   filename = "output_plots/hallmark_venn_all.png")
dev.off()
grid::grid.draw(vp)


# KEGG gene set enrichment ------------------------------------------------
# Full arguments
# kk2 <- gseKEGG(geneList = kegg_gene_list,
#                organism = kegg_org,
#                minGSSize = 3,
#                maxGSSize = 800,
#                pvalueCutoff = 0.05,
#                pAdjustMethod = "none",
#                keyType = "ncbi-geneid")
kegg_org <- "gfr"

# 1. fortis infected vs. uninfected
for_iu_kegg <- sort_genes(for_iu[for_iu$sex.linked=="-",], species = "fortis") %>%
  gseKEGG(., kegg_org, pAdjustMethod = "fdr")
for_iu_kegg  %>% ridgeplot(fill = "p.adjust")+ ggtitle("Ground: infected vs. uninfected")

# 2. cra infected vs uninfected
cra_iu_kegg <- sort_genes(cra_iu, species = "fortis") %>%
  gseKEGG(., kegg_org, pAdjustMethod = "fdr")
cra_iu_kegg  %>% ridgeplot(fill = "qvalue")+  ggtitle('Vegetarian: infected vs. uninfected')
#cra_iu_kegg  %>% gseaplot2(., geneSetID = 1:nrow(.))
#cra_iu_kegg@result %>% pull(Description)#$filter(Description %in% "Glutamine")

# 3. fortis infected vs. cra infected
cra_for_i_kegg <- sort_genes(cra_for_i, species = "fortis") %>%
  gseKEGG(., kegg_org, pAdjustMethod = "fdr")
cra_for_i_kegg %>% ridgeplot(fill = "p.adjust") + ggtitle("Infected: ground vs. vegetarian ")


# 4. fortis uninfected vs. cra uninfected
cra_for_u_kegg <- sort_genes(cra_for_u, species = "fortis") %>%
  gseKEGG(., kegg_org, pAdjustMethod = "fdr")
ridgeplot(cra_for_u_kegg, fill = "p.adjust") +  ggtitle('U fortis vs. U cra')
cra_for_u_gsea %>% gseaplot2(., geneSetID = 1:nrow(.))

# Now males vs. females
# fortis infected vs. uninfected males
for_miu_kegg <- sort_genes(for_miu, species = "fortis") %>%
  gseKEGG(., kegg_org, pAdjustMethod = "fdr")

#sort_genes(for_miu, species = "fortis") %>%
#  gseMKEGG(., kegg_org, pAdjustMethod = "none") %>% ridgeplot(., fill = "qvalue")

# fortis females infected vs uninfected
for_fiu_kegg <- sort_genes(for_fiu, species = "fortis") %>%
  gseKEGG(., kegg_org, pAdjustMethod = "fdr")
for_fiu_kegg %>% ridgeplot(fill = "p.adjust") + ggtitle("fortis female i vs u")
for_fiu_kegg %>% gseaplot2(., geneSetID = 1:nrow(.))

# cra males infected vs uninfected.
cra_miu_kegg <- sort_genes(cra_miu, species = "fortis") %>%
  gseKEGG(., kegg_org, pAdjustMethod = "fdr") #nothing significant
#cra_miu_fiu_kegg %>% ridgeplot(fill = "p.adjust") # nothing significant

# for vs cra infected
cra_for_mi_kegg <- sort_genes(for_cra_mi, species = "fortis") %>%
  gseKEGG(., kegg_org, pAdjustMethod = "fdr")

# for vs. cra uninfected
cra_for_mu_kegg <- sort_genes(for_cra_mu, species = "fortis") %>%
  gseKEGG(., kegg_org, pAdjustMethod = "fdr")

# kegg ridgeplots ---------------------------------------------------------

p8 <- for_iu_kegg %>% ridgeplot(fill = "p.adjust")+ ggtitle("Ground: infected vs. uninfected")+
  xlab(expression(log[2]*" fold change" ))
p9 <- cra_iu_kegg  %>% ridgeplot(fill = "qvalue")+  ggtitle('Vegetarian: infected vs. uninfected')+
  xlab(expression(log[2]*" fold change" ))
p10 <- cra_for_i_kegg %>% ridgeplot(fill = "p.adjust") + ggtitle("Infected: ground vs. vegetarian ")+
  xlab(expression(log[2]*" fold change" ))
p11 <- ridgeplot(cra_for_u_kegg, fill = "p.adjust") +  ggtitle("Uninfected: ground vs. vegetarian")+
  xlab(expression(log[2]*" fold change" ))

pdf("output_plots/kegg_ridgeplots.pdf", width = 12, height = 9)
plot_grid(p8, p9, p10, p11, labels = c('A', 'B', 'C', 'D'), label_size = 14)
dev.off()


# Now males vs. females
# fortis infected vs. uninfected males
p12 <- for_miu_kegg  %>% ridgeplot(fill = "p.adjust") + ggtitle("Ground: infected vs. uninfected males")+
  xlab(expression(log[2]*" fold change" ))
p13 <- cra_for_mi_kegg %>% ridgeplot(fill = "p.adjust") + ggtitle("Infected: male ground vs. male vegetarian")+
  xlab(expression(log[2]*" fold change" ))
p14 <- cra_for_mu_kegg %>% ridgeplot(fill = "p.adjust") + ggtitle("Uninfected: male ground vs. male vegetarian")+
  xlab(expression(log[2]*" fold change" ))

pdf("output_plots/kegg_ridgeplots_males.pdf", width = 12, height = 9)
plot_grid(p12,p13,p14, labels = c("A", "B", "C"), label_size = 14)
dev.off()

# Kegg venn diagram. Not sure this is necessary because there are such few categories
# involved here.
vp_kegg <- venn.diagram(list (
  "Uninf. vs. uninf.\nN = 3" = cra_for_u_kegg@result$ID,
  "Inf. vs inf.\nN = 8" = cra_for_i_kegg@result$ID,
  "Vegetarian\nN = 9" = cra_iu_kegg@result$ID,
  "Ground\nN = 14" = for_iu_kegg@result$ID),
  fill = 2:5, alpha = 0.3, cat.cex = 1.5, cex =1.4,
  #print.mode = "percent", sigdigs = 2,
  filename = "output_plots/kegg_venn_all.png")
dev.off()
grid::grid.draw(vp_kegg)


#
# Cluster Profiler other functions ---------------------------------------------------

# http://yulab-smu.top/clusterProfiler-book/chapter2.html#gene-set-enrichment-analysis

# 1. Universal enrichment analysis:
# hypergeometric tests through enricher() and GSEA() for gene set enrichment
# Over representation analysis just starts with a vector of gene IDS i.e.
# the output of DESeq2
#
# Wikipathways analysis:
# Download gmt file from data.wikipathways.org
# http://data.wikipathways.org/current/gmt/wikipathways-20210210-gmt-Gallus_gallus.gmt

# Data (male for i vs u)
ids <- queryMany(genes$GeneID, scope = "symbol", fields = "entrezgene", species = "48883") %>%
  as.data.frame %>% dplyr::select(query, entrezgene) %>% filter(!is.na(entrezgene)) %>%
  rename(GeneID = query)

ids_gallus <- queryMany(genes$GeneID, scope = "symbol", fields = "entrezgene", species = "9031") %>%
  as.data.frame %>% dplyr::select(query, entrezgene) %>% filter(!is.na(entrezgene)) %>%
  rename(GeneID = query)


head(genes)
# Pull out the list of NCBI gene IDs where abs(log2 fold change) > 2
wiki_gene_list <- filter(genes, GeneID %in% ids_gallus$GeneID ) %>%
  merge(ids_gallus) %>%
  filter(abs(log2.FC.) > 2) %>%
  pull(entrezgene)

# pull and clean wiki pathways gmt file for chicken
wp2gene <- read.gmt("input_data/wikipathways-20210210-gmt-Gallus_gallus.gmt")
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
ewp <- enricher(wiki_gene_list, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
head(ewp) # nothing comes out significant

# Nothing significant here
GSEA_chickn <- filter(genes, GeneID %in% ids_gallus$GeneID ) %>%
  merge(ids_gallus) %>%
  select(log2.FC., "entrezgene") %>%
  arrange(desc(log2.FC.)) %>% pull(log2.FC., name = entrezgene) %>%
  GSEA(.,  TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose = T)


# Section 5: Gene ontology analysis
# groupGO(), enrichGO() and gseGO()

# groupGO: classification of genes for a specific level of GO hierarchy
GOlist <- filter(genes, GeneID %in% ids_gallus$GeneID ) %>%
  merge(ids_gallus) %>%
  filter(abs(log2.FC.) > 2) %>% pull(entrezgene)
groupGO(gene = GOlist,
        OrgDb = org.Gg.eg.db,
        ont = "BP",
        readable = T,
        level = 3 ) -> ggo
ggo[1:30]

# GO over-representation test
ego <- enrichGO(gene = GOlist,
                universe = names(ids_gallus$entrezgene),
                OrgDb = org.Gg.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)

#no enriched terms found with adjusted P
ego

# GO Gene Set Enrichment Analysis
filter(genes, GeneID %in% ids_gallus$GeneID ) %>%
  merge(ids_gallus) %>%
  select(log2.FC., "entrezgene") %>%
  arrange(desc(log2.FC.)) %>% pull(log2.FC., name = entrezgene) %>%

  gseGO(geneList = .,
        OrgDb =  org.Gg.eg.db,
        ont = "BP",
        #nperm = 1000,
        minGSSize = 100,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        verbose = TRUE) -> ego3 # no enriched terms found





# KEGG --------------------------------------------------------------------

# KEGG analysis
# Another way to visualize pathways

# Translating key types, can use clusterProfiler::keytypes(), but relies
# on chicken database. mygene::quereyMany() will work with fortis directly.


# My gene style
ids <- queryMany(genes$GeneID, scope = "symbol", fields = "entrezgene", species = "48883") %>%
  as.data.frame %>% dplyr::select(query, entrezgene) %>% filter(!is.na(entrezgene)) %>%
  rename(GeneID = query)

kegg_gene_list <- filter(genes, GeneID %in% ids$GeneID ) %>%
  merge(ids) %>%
  select(log2.FC., "entrezgene") %>%
  arrange(desc(log2.FC.)) %>% pull(log2.FC., name = entrezgene)
head(kegg_gene_list)
#"gfr" = ground finch
kegg_org <- "gfr"


# KEGG over-representation test
# no overrepresented genes found
filter(genes, GeneID %in% ids$GeneID ) %>%
  merge(ids) %>%
  filter(log2.FC. > 2) %>%
  pull("entrezgene") %>%
  enrichKEGG(gene = .,
             organism = kegg_org,
             pvalueCutoff = 0.05)


# GSE KEGG
#
kegg_gene_list <- filter(genes, GeneID %in% ids$GeneID ) %>%
  merge(ids) %>%
  select(log2.FC., "entrezgene") %>%
  arrange(desc(log2.FC.)) %>% pull(log2.FC., name = entrezgene)
kk2 <- gseKEGG(geneList = kegg_gene_list,
               organism = kegg_org,
               # nperm = 1000, #apparently deprecated
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType = "ncbi-geneid")
kk2
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
pairwise_termsim(kk2) %>% emapplot
ridgeplot(kk2)

gseMKEGG(kegg_gene_list,
         organism = kegg_org,
         pAdjustMethod = "none")

# Gene Set Enrichment- FOR I vs. U ----------------------------------------
# Based on https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

# Next: read in data, e.g. fortis infected vs. uninfected
# Data frame is all 16000 genes, 1 row per gene
genes <- read.csv("input_data/for_m_i_u_clusterprofiler.csv") # just males

original_gene_list <- genes$log2.FC.  # Focus on log 2 fold change, i.e. how different expression was
names(original_gene_list) <- genes$GeneID # Name the genes
gene_list <- na.omit(original_gene_list) # Remove NAs
gene_list <- sort(gene_list, decreasing = TRUE) # sort from biggest to smallest

gse <- gseGO(geneList = gene_list,
             ont = "BP",
             keyType = "SYMBOL",
             #nPerm = 100, #apparently deprecated
             minGSSize = 3, # min number of genes in gene set
             maxGSSize = 800, # max number of genes in gene set
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = get("org.Gg.eg.db"), #avoids conflicts in species() fun
             pAdjustMethod = "none",
             by = "fgsea" # arbitrary choice between that and DOSE
)
gse$Description
dotplot(gse, showCategory=10, split=".sign") +
  facet_grid(.~.sign) #wtf
pairwise_termsim(gse) %>% emapplot() #create a network

#gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

ridgeplot(gse) + labs(x = "enrichment distribution") # make a ridgeplot






# ClusterProfiler_reactome ------------------------------------------------

# Download immune gene sets for chicken
reactome <- msigdbr(species = "Gallus gallus", category = "C2", subcategory = "CP:REACTOME")
reactome <- reactome %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
head(reactome)

for_iu_gsea <- sort_genes(for_iu[for_iu$sex.linked=="-",], species = "gallus") %>%
  GSEA(., TERM2GENE = reactome, pAdjustMethod = "none")
for_iu_gsea %>% dim # 14 sig dif
for_iu_gsea %>% ridgeplot()#label_format = function(x) uncap(x))+
  ggtitle("Ground: Infected vs uninfected") +
  xlab(expression(log[2]*" fold change" ))

# 2. fortis infected vs. cra infected
cra_for_i_gsea <- sort_genes(cra_for_i[cra_for_i$sex.linked=="-",], species = "gallus") %>%
  GSEA(., TERM2GENE = msigdbr_2gene, pAdjustMethod = "fdr")
ridgeplot(cra_for_i_gsea, label_format = function(x) uncap(x)) +  ggtitle('ground vs. veggie, all infected')
cra_for_i_gsea@result %>% dim # 9 sig diff
cra_for_i_gsea %>% gseaplot2(., geneSetID = 1:nrow(.))

# 3. cra infected vs uninfected
cra_iu_gsea <- sort_genes(cra_iu[cra_iu$sex.linked=="-",], species = "gallus") %>%
  GSEA(., TERM2GENE = msigdbr_2gene, pAdjustMethod = 'fdr')
cra_iu_gsea %>% ridgeplot(label_format = function(x) uncap(x))+  ggtitle('Veggie: infected vs. uninfected')
#cra_iu_gsea %>% gseaplot2(., geneSetID = 1:nrow(.))



# 4. fortis uninfected vs. cra uninfected
cra_for_u_gsea <- sort_genes(cra_for_u[cra_for_u$sex.linked=="-",], species = "gallus") %>%
  GSEA(., TERM2GENE = msigdbr_2gene, pAdjustMethod = 'fdr')
