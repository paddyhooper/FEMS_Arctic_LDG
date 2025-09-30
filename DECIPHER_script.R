##Script for creating phylogenetic trees in DECIPHER
#Author: Patrick Hooper
#Date Created: 23/02/2022

#This script requires installation of the r_env conda environment 
#This script can be used for both the 16S and 18S rRNA phyloseq packages

#install these on conda
#   conda install -c bioconda r-phangorn 
#   conda install -c bioconda bioconductor-phyloseq
#   conda install -c bioconda bioconductor-decipher

#Load libraries
library(phyloseq); packageVersion("phyloseq")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")

#IMPORTANT: You will need to direct the script to your desired phyloseq object

print("ready to begin!")
#step 1
ps <- readRDS(file = "~/filt_ps_18S")
ps

filterMetazoa_NA = c("Metazoa","NA") #change to your specific Supergroup(s)

# Remove filtered Supergroups
ps = subset_taxa(ps, !Division %in% filterMetazoa_NA)
ps

print("beginning tree assembly")

# Phylogenetic tree ####
seqs <- DNAStringSet(refseq(ps))
seqs
alignment <- AlignSeqs(seqs, anchor=NA)
print("tree alignment complete")

# Phangorn package ####
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

dm <- dist.ml(phang.align)

saveRDS(dm, file = "~/dm")
print("tree data matrix saved")

treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

saveRDS(fitGTR, file = ~/fitGTR)
ps.merged <- merge_phyloseq(ps, phy_tree(fitGTR$tree))

print("tree added to phyloseq object")

# Save the final test phyloseq object with tree included
saveRDS(ps.merged_test, "~/ps.merged")

print("final phyloseq object saved")



