#Script for running tip_glom
#Author: Patrick Hooper
#Date Created: 02/03/2022

library(phyloseq); packageVersion("phyloseq")
print("ready to begin!")

ps <- readRDS(file = "~/ps_merged_16S_2_4_22")
ps

h2 = 0.02
print("glomming 0.02")

ps_merged_16S_glom_0.02 = tip_glom(ps, h = h2)
ps_merged_16S_glom_0.02
saveRDS(ps_merged_16S_glom_0.02, file = "~/canada/18S/ps_merged_16S_glom_0.02_4_4")

print("glom 0.02 finished")

h1 = 0.01
print("glomming 0.01")

ps_merged_16S_glom_0.01 = tip_glom(ps, h = h1)
ps_merged_16S_glom_0.01
saveRDS(ps_merged_16S_glom_0.01, file = "~/canada/18S/ps_merged_16S_glom_0.01_4_4")

print("glom 0.01 finished")

h3 = 0.2
print("glomming 0.2")

ps_merged_16S_glom_0.2 = tip_glom(ps, h = h3)
ps_merged_16S_glom_0.2
saveRDS(ps_merged_16S_glom_0.2, file = "~/canada/18S/ps_merged_16S_glom_0.2_4_4")

print("Fin.")

