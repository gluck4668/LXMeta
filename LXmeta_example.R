
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

devtools::load_all()

file_data="meta_data.xlsx"

species="mouse"  # The species should be human, mouse, or rat.


LXmeta(file_data,species)



#---If MetaboAnalystR_3.2.0 can not be installed properly, please use the meature below:

# step1: download the package from : https://www.dropbox.com/s/pp9vziji96k5z5k/MetaboAnalystR_3.2.0.tar.gz
# step2: install.packages("D:/MetaboAnalystR_3.2.0.tar.gz", repos = NULL, type = "source")

