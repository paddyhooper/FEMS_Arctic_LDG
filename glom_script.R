#Script for phylogenetic agglomeration of 16S and 18S ASVs
#Author: Patrick M. Hooper
#Date Created: 02/03/2022

#This script is run in command line and relies on having an activated r-env

library(phyloseq); packageVersion("phyloseq")
print("ready to begin!")

#Direct this to your 16S and 18S phyloseq objects
ps <- readRDS(file = "~/ps_merged")
ps

h1 = 0.01
print("glomming 0.01")

ps_merged_16S_glom_0.01 = tip_glom(ps, h = h1)
ps_merged_16S_glom_0.01
saveRDS(ps_merged_16S_glom_0.01, file = "~/canada/18S/ps_merged_16S_glom_0.01")

print("glom 0.01 finished")

#End of Script
