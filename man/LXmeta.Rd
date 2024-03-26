
\name{LXmeta}
\alias{LXmeta}
\title{LXmeta: metabolite set enrichment and pathways analysis}
\usage{
LXmeta(file_data,species)
}
\description{
This R package can be used to analyze the metabolite enrichment and the related metabolism pathways.
It cited the R packages publised by "https://www.metaboanalyst.ca/MetaboAnalyst".
}
\examples{

if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXmeta")

library(LXmeta)

??LXmeta

#-----------------------
data(meta_data_example)

#-----------------------

rm(list=ls())

file_data="meta_data.xlsx"

species="rat"  # The species should be human, mouse, or rat.


LXmeta(file_data,species)

}
