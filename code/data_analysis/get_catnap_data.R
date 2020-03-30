#!/usr/local/bin/Rscript

## load required functions and packages
library("RCurl")

## -------------------------------------------------------------------------
## get data from Magaret et al. (2019)
## -------------------------------------------------------------------------
## load data set 1
data1 <- read.csv(text = getURL("https://raw.githubusercontent.com/benkeser/vrc01/master/data/data1.csv"), header = TRUE)
# load data set 2
data2 <- read.csv(text = getURL("https://raw.githubusercontent.com/benkeser/vrc01/master/data/data2.csv"), header = TRUE)

## merge them together
data1$dataset <- 1
data2$dataset <- 2
dataset <- rbind(data1, data2)
if (!dir.exists("../code/data/")) {
    dir.create("../code/data/", recursive = TRUE)
}
saveRDS(dataset, "data/analysis_data.rds")
