
library(tidyverse)
library(dplyr)



# Data --------------------------------------------------------------------
metadata <- read.csv("input_data/metadata.csv")
bracken <- read.table("bracken/bracken.raw.viral.0.05.export.text", sep = "\t",
                      col.names = c("sample2", "family", "taxonid", "level",
                                    "kraken_assigned", "added_reads",
                                    "new_est_reads", "fraction_total_reads", "analysis"))
with(metadata[metadata$sex=="M",], table(species, infection.status))
# Output from kraken classifying using protein database
# Taxon id for pox family: Chordopoxviriniae 10241
# Taxon id for avipoxvirus genus: 10260
kraken_protein <- read.table("bracken/kraken.protein.0.00.export.txt", sep = "\t",
                             col.names = c("sample2", "percent.reads",
                                           "number.frag.group", "number.frag",
                                           "level", "taxonid", "taxon"))
kraken_protein  <- kraken_protein %>% filter(taxonid == "10241") %>% merge(., metadata)

kraken_protein %>% filter(taxonid == "10260") %>% merge(., metadata) %>%
aggregate(number.frag.group ~ species + infection.status, ., mean)
fivenum(kraken_protein$number.frag.group)

glm.nb(number.frag.group/unmapped_reads ~ species + infection.status, data = kraken_protein) %>%
        summary()
var(kraken_protein$number.frag.group)
aggregate(number.frag.group ~ + infection.status + species , data = kraken_protein, mean)

kraken_protein_5 <- read.table("bracken/kraken.protein.0.05.export.txt", sep = "\t",
                             col.names = c("sample2", "percent.reads",
                                           "number.frag.group", "number.frag",
                                           "level", "taxonid", "taxon"))

kraken_protein_5  <- kraken_protein_5 %>% filter(taxonid == "10241") %>% merge(., metadata, all.y = T) #%>%
        aggregate(number.frag.group ~ species + infection.status, ., mean)

glm(number.frag.group ~ species + infection.status, data = kraken_protein_5, family = "poisson") %>%
        summary()


for_miu <- read.csv("input_data/for_m_i_u_clusterprofiler.csv") # fortis infected vs not
for_iu <- read.csv("input_data/for_i_u_clusterprofiler.csv")
counts <- read.csv("input_data/normalized_counts_dseq_all.csv")
raw_counts <- read.table("input_data/2008D-Gfortis_rawCounts.txt", row.names = NULL)
head(raw_counts)




# Summary stats -----------------------------------------------------------
# Some numbers about sequencing results and alignment etc.

metadata$perc.align <- metadata$mapped.reads/metadata$total_reads
fivenum(metadata$perc.align)

summary(lm(mapped.reads ~ species, metadata))

# Bracken -----------------------------------------------------------------

head(bracken)
bracken <- merge(bracken, metadata)
bracken <- bracken %>% filter(family == "Poxviridae")
boxplot(fraction_total_reads ~ species + infection.status, bracken,
        ylab = "Fraction of total unmapped reads",
        names = rep(c("Veggie", "Med. ground"), 2), xlab = ""
        )
mtext("Infected", side = 1, line = 3, adj = 0.25, cex = 1.2)
mtext("Uninfected", side = 1, line = 3, adj = 0.75, cex = 1.2)

lm(fraction_total_reads/unmapped_reads ~ species + infection.status, bracken) %>% summary

for_iu %>% filter(DIFF_EXP==1) %>% arrange(desc(avg.IFOR))

# special gene DDX58 av. expression infected 1860      449,
# another gene DHX58
#
ddx2 <-filter(raw_counts, row.names == "DHX58")[,-1] %>% t %>% as.data.frame
ddx2 <- cbind(ddx2, str_split(rownames(ddx2), pattern = "[.]", simplify = T ))[,c(1:2)]
colnames(ddx2) <- c("raw_reads", "sample2")
ddx2 <- merge(ddx2, bracken)

ddx <- filter(counts, X == "DDX58")[,-1] %>% t %>% as.data.frame
ddx <- cbind(ddx, str_split(rownames(ddx), pattern ="[.]", simplify = T))[,c(1:2)]
colnames(ddx) <- c("normalized_reads", "sample2")
ddx <- merge(ddx, bracken)

plot(fraction_total_reads ~ normalized_reads, ddx,
     col =  c("#117733", "#882255", "#44AA99", "#CC6677")[as.factor(group)], pch = 19,
     ylab = "Fraction of total reads assigned to pox", xlab = "Normalized DDX58 reads", cex = 1.2)
legend("topright", legend = c("I-veg", "I-ground", "U-veg", "U-ground"), pch = 19,
       col =c("#117733", "#882255", "#44AA99", "#CC6677"))
lm(fraction_total_reads ~ normalized_reads * species + , ddx) %>% summary
sessionInfo()
library(clusterProfiler)
