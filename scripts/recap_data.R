library(dplyr)
library(lme4)


# Data
hendry <- read.csv("banding_data/hendry_bands.csv")
birgit <- read.csv("banding_data/birgit_bands.csv")
sabrina <- read.csv("banding_data/sabrina_bands.csv")

head(hendry)
head(birgit)
head(sabrina)


# analysis

infected <- filter(sabrina, Species %in% c("CRA", "FOR") & POX %in% c("I", "I-R")) %>% droplevels()

infected$band
filter(hendry, band %in% infected$band)
filter(birgit, band %in% infected$band & Year..four.digitals. %in% c("2016","2017", "2018", "2019") )
