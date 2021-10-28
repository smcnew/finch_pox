# Script to analyze slide microscopy data

library(MASS)
library(dplyr)



# Data --------------------------------------------------------------------


slides <- read.csv("input_data/leukocyte_results.csv") %>%
  filter(Read_unread == "read") %>% filter(pox != "R") %>%
  mutate(hl_ratio = heterophil/lymphocyte) %>%
  mutate(clade = case_when(species == "FOR" ~ "ground",
                           species == "FUL" ~ "ground",
                           species == "CRA" ~ "veggie"))
slides$hl_ratio[slides$hl_ratio == "Inf"] <- NA

# Summary stats
filter(slides, species != "FUL") %>% dplyr::select(total.leukocyte, total.erythrocytes) %>%
apply(2, median) # median counts of leukos and erythrocytes
  # filter(., total.leukocyte < 50) %>% pull(total.erythrocytes)

fivenum(slides$total.leukocyte)
table(slides$species, slides$pox)

# Just veggie and fortis --------------------------------------------------
# filter out fuliginosa
slides2 <- filter(slides, species != "FUL") %>% droplevels()
fivenum(slides2$total.leukocyte)



# Plots -------------------------------------------------------------------

colors <- c("#009988", "#EE7733")
pdf("output_plots/leukocytes.pdf")
par(mfrow=c(2,2), mar = c(5,4,2,2))
boxplot(leukocyte.per.1000.erythrocyte ~ species + pox, data=  slides2,
        ylab = "Leukocytes/1000 erythrocytes",
        names = rep(c("Veg", "Ground"), 2),
        xlab = "", col = rep(colors, 2))
mtext("Infected", side = 1, line = 3, adj = 0.2, cex = 1.2)
mtext("Uninfected", side = 1, line = 3, adj = 0.8, cex = 1.2)
stripchart(leukocyte.per.1000.erythrocyte ~ species + pox, data=  slides2,
           add = T, vertical = T, pch = 21, bg = "grey30", cex = 0.7)
mtext("A", side = 3, line =0, adj = - .23, cex = 1.4)

boxplot(lymphocyte.per.1000.erythrocyte ~ species + pox, data=  slides2,
        ylab = "Lymphocytes/1000 erythrocytes",
        names = rep(c("Veg", "Ground"), 2),
        xlab = "", col = rep(colors, 2))
mtext("Infected", side = 1, line = 3, adj = 0.2, cex = 1.2)
mtext("Uninfected", side = 1, line = 3, adj = 0.8, cex = 1.2)
stripchart(lymphocyte.per.1000.erythrocyte ~ species + pox, data=  slides2,
           add = T, vertical = T, pch = 21, bg = "grey30", cex = 0.7)
mtext("B", side = 3, line = 0, adj = - .23, cex = 1.4)

boxplot(monocyte.per.1000.erythrocyte ~ species + pox, data=  slides2,
        ylab = "Monocytes/1000 erythrocytes",
        names = rep(c("Veg", "Ground"), 2),
        xlab = "", col = rep(colors, 2))
mtext("Infected", side = 1, line = 3, adj = 0.2, cex = 1.2)
mtext("Uninfected", side = 1, line = 3, adj = 0.8, cex = 1.2)
stripchart(monocyte.per.1000.erythrocyte ~ species + pox, data=  slides2,
           add = T, vertical = T, pch = 21, bg = "grey30", cex = 0.7)
mtext("C", side = 3, line = 0, adj = - .23, cex = 1.4)

boxplot(log(hl_ratio+1) ~ species + pox, data=  slides2,
        ylab = "log(Heterophil:lymphocyte ratio) + 1 ",
        names = rep(c("Veg", "Ground"), 2),
        xlab = "", col = rep(colors, 2))
stripchart( log(hl_ratio+1) ~ species + pox, data=  slides2,
            add = T, vertical = T, pch = 21, bg = "grey30", cex = 0.7)
mtext("Infected", side = 1, line = 3, adj = 0.2, cex = 1.2)
mtext("Uninfected", side = 1, line = 3, adj = 0.8, cex = 1.2)
mtext("D", side = 3, line = 0, adj = - .23, cex = 1.4)

dev.off()

# Stats -------------------------------------------------------------------


summary(glm(leukocyte.per.1000.erythrocyte ~ species + pox, family = "quasipoisson", data=  slides2))
summary(glm(lymphocyte.per.1000.erythrocyte ~ species, family = "quasipoisson", data=  slides2))
summary(glm(monocyte.per.1000.erythrocyte ~ species + pox, family = "quasipoisson", data=  slides2))
summary(glm(hl_ratio ~ species + pox , family = "quasipoisson", data=  slides2))

# All species -------------------------------------------------------------


# All Leukos
colors <- c("#009988", "#EE7733", "#CC3311")
boxplot(leukocyte.per.1000.erythrocyte ~ species + pox, data=  slides,
        ylab = "Leukocytes per 1000 erythrocytes",
        names = rep(c("Veggie", "Med. ground", "Sm. ground"), 2),
        xlab = "", col = rep(colors, 2))
mtext("Infected", side = 1, line = 3, adj = 0.25, cex = 1.2)
mtext("Uninfected", side = 1, line = 3, adj = 0.75, cex = 1.2)
stripchart(leukocyte.per.1000.erythrocyte ~ species + pox, data=  slides,
           add = T, vertical = T, pch = 21, bg = "grey30", cex = 0.7)

summary(lm(leukocyte.per.1000.erythrocyte ~ species, data=  slides))

boxplot(lymphocyte.per.1000.erythrocyte ~ species + pox, data=  slides,
        ylab = "Lymphocytes per 1000 erythrocytes",
        names = rep(c("Veggie", "Med. ground", "Sm. ground"), 2),
        xlab = "", col = rep(colors, 2))
stripchart(lymphocyte.per.1000.erythrocyte ~ species + pox, data=  slides,
           add = T, vertical = T, pch = 21, bg = "grey30", cex = 0.7)
mtext("Infected", side = 1, line = 3, adj = 0.25, cex = 1.2)
mtext("Uninfected", side = 1, line = 3, adj = 0.75, cex = 1.2)

summary(lm(lymphocyte.per.1000.erythrocyte ~ clade * pox, data=  slides))


boxplot(heterophil.per.1000.erythrocyte ~ species + pox, data=  slides,
        ylab = "Heterophils per 1000 erythrocytes",
        names = rep(c("Veggie", "Med. ground", "Sm. ground"), 2),
        xlab = "", col = rep(colors, 2))
stripchart(lymphocyte.per.1000.erythrocyte ~ species + pox, data=  slides,
           add = T, vertical = T, pch = 21, bg = "grey30", cex = 0.7)
mtext("Infected", side = 1, line = 3, adj = 0.25, cex = 1.2)
mtext("Uninfected", side = 1, line = 3, adj = 0.75, cex = 1.2)

summary(lm(lymphocyte.per.1000.erythrocyte ~ species + pox, data=  slides))


boxplot(log(hl_ratio+1) ~ species + pox, data=  slides,
        ylab = "log(Heterophil:lymphocyte ratio) + 1 ",
        names = rep(c("Veggie", "Med. ground", "Sm. ground"), 2),
        xlab = "", col = rep(colors, 2))
stripchart( log(hl_ratio+1) ~ species + pox, data=  slides,
           add = T, vertical = T, pch = 21, bg = "grey30", cex = 0.7)
mtext("Infected", side = 1, line = 3, adj = 0.25, cex = 1.2)
mtext("Uninfected", side = 1, line = 3, adj = 0.75, cex = 1.2)

summary(lm(log(hl_ratio + 1) ~ species + pox + sex, data=  slides))
summary(glm.nb(hl_ratio+1 ~ species + pox , data=  slides))

