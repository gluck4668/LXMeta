
setwd("D:/R-lin study/R packages/LXmeta")
library(openxlsx)

meta_data_example <- read.xlsx("meta_data.xlsx")

usethis::use_data(meta_data_example,overwrite = T)

rm(list=ls())

data(meta_data_example)

