#Assembly of V4 16S sequence data using dada2 ASV assembly 
#Author: Patrick Hooper
##Created: 01.04.22

#Along the way I used saveRDS to save R products, this can save a lot of time in the future as you can just ##Reload previously run command outputs using the script in my GitHub page: "Reloading your DADA2 analyses"

#OPTIONAL: Add the MOTHUR SOP mock community to your 16S reads to test the accuracy of your assignment
#Available: https://benjjneb.github.io/dada2/tutorial.html, section 'Getting Ready' from hyper link 'example data'
#Copy across the F1 and R1 mock community '.fa' files and the reference file 'HMP_Mock.v35.fasta' to the folder containing your sequencing files

#INSTALLATION####
install.packages("devtools")
devtools::install_github("benjjneb/dada2", ref = "v1.16") # change the ref argument to get other versions
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
devtools::install_github("KlausVigo/phangorn")

#LOAD PACKAGES####
library(devtools); packageVersion("devtools")
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(tidyverse); packageVersion("tidyverse")
library(ggplot2); packageVersion("ggplot2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(microbiome); packageVersion("microbiome")
library(knitr); packageVersion("knitr")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
library(gridExtra);packageVersion("gridExtra")
theme_set(theme_minimal()) # set ggplot2 theme to minimal

#SET YOUR WORKING DIRECTORY
setwd("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2")

#YOU'RE READY TO BEGIN WITH DADA2!
#My sequencing files are from two separate sequencing runs, DADA2 requires you to assemble your two ASV tables separately and merge them after assembly
#I have two sequencing runs from 2016 and 2018, i'll start with the 2016 sequences

#PART 1: IDENTIFY PATH TO YOUR SEQUENCING FILES
path <- "C:/Users/pmh36/OneDrive - Natural History Museum/Data/Canada_Datasets/Qiime_16S_2016/2016_16S_Sequences/"
list.files(path) # check your files are all present

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs_16S_2016 <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs_16S_2016 <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs_16S_2016), "_"), `[`, 1)
sample.names # check you've got all your samples (and your mock community if using)!

saveRDS(fnFs_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/fnFs_16S_2016")
saveRDS(fnRs_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/fnRs_16S_2016")

#SECTION 1:
#Inspect read quality profiles for FWD reads for first 4 sequences
FWD_Quality_Profile <- plotQualityProfile(fnFs_16S_2016[1:4])
#Inspect read quality profiles for REV reads for first 4 sequences
REV_Quality_Profile <- plotQualityProfile(fnRs_16S_2016[1:4])
FWD_Quality_Profile
REV_Quality_Profile

saveRDS(FWD_Quality_Profile, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/FWD_Quality_Profile")
saveRDS(REV_Quality_Profile, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/REV_Quality_Profile")

#Grey scale heat map is the frequency of each quality score at each base position
#Mean quality score = green line
#Quartiles of quality score = orange line
#Red line = scaled proportion of reads extend to at least that position, not useful for illumina as all reads same length

#SECTION 2: Filter and Trim
# Place filtered files in filtered/ subdirectory
filtFs_16S_2016 <-
  file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_16S_2016 <-
  file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs_16S_2016) <- sample.names
names(filtRs_16S_2016) <- sample.names

# In most microbial species,the 16S fourth hypervariable (V4) region is approximately 254 bp, and only deviates from this length by a few basepairs. 
#Worth considering how best to trim to produce the best results
out_16S_2016 <- filterAndTrim(
  fnFs_16S_2016,
  filtFs_16S_2016,
  fnRs_16S_2016,
  filtRs_16S_2016,
  truncLen = c(220, 150), #truncate reads after this many base pairs
  maxN = 0,  #after truncation, reads with more than maxNs will be discarded
  maxEE = c(2, 2), #maximum number of Expected Errors in a read. If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5))
  truncQ = 2, #truncates reads at the first instance of a quality score less than or equal to trunq
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE #not possible on windows
) # On Windows set multithread = FALSE
out_16S_2016

saveRDS(out_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/out_16S_2016")

#SECTION 3: Learn the error rates
errF_16S_2016 <- learnErrors(filtFs_16S_2016, multithread = TRUE)
errR_16S_2016 <- learnErrors(filtRs_16S_2016, multithread = TRUE)

saveRDS(errF_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/errF_16S_2016")
saveRDS(errR_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/errR_16S_2016")
#The error rates for each possible nucleotide transition are shown. Points = error rates for each quality score, black line = estimated error rates from the machine-learning algorithm. Red lines shows expected error rate under nominal Q-score. If points fit well to the black line and error rate drops with increased quality score you can assume a good error rate.
Plot_Error_FWD <- plotErrors(errF_16S_2016, nominalQ = TRUE)#FWD read errors
Plot_Error_REV <- plotErrors(errR_16S_2016, nominalQ = TRUE) #REV read errors
Plot_Error_FWD
Plot_Error_REV

saveRDS(Plot_Error_FWD, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/Plot_Error_FWD")
saveRDS(Plot_Error_REV, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/Plot_Error_REV")

#SECTION 4: Sample inference
#Apply the core sample inference algorithm to the filtered and trimmed sequence data
#Info: https://www.nature.com/articles/nmeth.3869#methods
dadaFs_16S_2016 <- dada(filtFs_16S_2016, err = errF_16S_2016, multithread = FALSE)
dadaRs_16S_2016 <- dada(filtRs_16S_2016, err = errR_16S_2016, multithread = FALSE)

saveRDS(dadaFs_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/dadaFs_16S_2016")
saveRDS(dadaRs_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/dadaRs_16S_2016")

#Inspect the first line
dadaFs_16S_2016[[1]]

#4.2. Merge your filtered sequence files with your sample inference
mergers_16S_2016 <- mergePairs(dadaFs_16S_2016, filtFs_16S_2016, dadaRs_16S_2016, filtRs_16S_2016, verbose = TRUE)

# Inspect the merger data.frame from the first sample
head(mergers_16S_2016[[1]])

saveRDS(mergers_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/mergers_16S_2016")

#CONSTRUCT SEQUENCE TABLE
seqtab_16S_2016 <- makeSequenceTable(mergers_16S_2016)
dim(seqtab_16S_2016)
saveRDS(seqtab_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/seqtab_16S_2016")

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_16S_2016)))

# Remove all sequences outside the length range of 250-256.
seqtab_16S_2016_trim <- seqtab_16S_2016[,nchar(colnames(seqtab_16S_2016)) %in% seq(250,256)]
table(nchar(getSequences(seqtab_16S_2016_trim)))

saveRDS(seqtab_16S_2016_trim, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/seqtab_16S_2016_filtered")

#REMOVE CHIMERAS
seqtab_nochim_16S_2016 <- removeBimeraDenovo(seqtab_16S_2016_trim, method = "consensus", multithread = FALSE, verbose = TRUE)
dim(seqtab_nochim_16S_2016)

saveRDS(seqtab_nochim_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/seqtab_nochim_16S_2016")

# Inspect distribution of sequence lengths
seq_length_nochim_16S_2016 <- table(nchar(getSequences(seqtab_nochim_16S_2016)))
seq_length_nochim_16S_2016

write.csv(seq_length_nochim_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/seq_length_distribution")

ncol(seqtab_16S_2016)
ncol(seqtab_nochim_16S_2016)
sum(seqtab_nochim_16S_2016) / sum(seqtab_16S_2016) * 100 # % of non-chimeras

#SECTION 5: Track reads through the pipeline
getN <- function(x)
  sum(getUniques(x))
getN

track_16S_2016 <-
  cbind(
    out_16S_2016,
    sapply(dadaFs_16S_2016, getN),
    sapply(dadaRs_16S_2016, getN),
    sapply(mergers_16S_2016, getN),
    rowSums(seqtab_nochim_16S_2016)
  )
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track_16S_2016) <-
  c("input",
    "filtered",
    "denoisedF",
    "denoisedR",
    "merged",
    "nonchim")
rownames(track_16S_2016) <- sample.names
print(track_16S_2016)

write.table(track_16S_2016, "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/track_16S_2016.txt")
saveRDS(track_16S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/track_16S_2016")

#Outside of filtering, there should no step in which a majority of reads are lost
#If a majority of reads failed to merge, you may need to revisit the truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon. #If a majority of reads were removed as chimeric, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification

#OPTIONAL PART 2: IDENTIFY PATH TO YOUR SECOND FOLDER OF SEQUENCING FILES
path <- "C:/Users/pmh36/OneDrive - Natural History Museum/Data/Canada_Datasets/Qiime_16S_2018/2018_16S_Sequences/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs_16S_2018 <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs_16S_2018 <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

saveRDS(fnFs_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/fnFs_16S_2018")
saveRDS(fnRs_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/fnRs_16S_2018")
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names_16S_2018 <- sapply(strsplit(basename(fnFs_16S_2018), "_"), `[`, 1)
sample.names_16S_2018

#SECTION 1: Inspect read quality profiles for FWD reads for first 4 sequences
FWD_Quality_Profile_2018 <- plotQualityProfile(fnFs_16S_2018[1:4])
#Inspect read quality profiles for REV reads for first 4 sequences
REV_Quality_Profile_2018 <- plotQualityProfile(fnRs_16S_2018[1:4])

FWD_Quality_Profile_2018
REV_Quality_Profile_2018

saveRDS(FWD_Quality_Profile_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/FWD_Quality_Profile")
saveRDS(REV_Quality_Profile_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/REV_Quality_Profile")

#Grey scale heat map is the frequency of each quality score at each base position
#Mean quality score = green line
#Quartiles of quality score = orange line
#Red line = scaled proportion of reads extend to at least that posotion, not useful for illumina as all reads same length

#SECTION 2: Filter and Trim
# Place filtered files in filtered/ subdirectory
filtFs_16S_2018 <- file.path(path, "filtered", paste0(sample.names_16S_2018, "_F_filt.fastq.gz"))
filtRs_16S_2018 < file.path(path, "filtered", paste0(sample.names_16S_2018, "_R_filt.fastq.gz"))

names(filtFs_16S_2018) <- sample.names_16S_2018
names(filtRs_16S_2018) <- sample.names_16S_2018
names(filtRs_16S_2018)
names(filtFs_16S_2018)

out_16S_2018 <- filterAndTrim(
  fnFs_16S_2018,
  filtFs_16S_2018,
  fnRs_16S_2018,
  filtRs_16S_2018,
  truncLen = c(230, 150),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE
) # On Windows set multithread=FALSE
out_16S_2018

saveRDS(out_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/out_16S_2018")


#SECTION 3: Learn the error rates
errF_16S_2018 <- learnErrors(filtFs_16S_2018, multithread = TRUE)
errR_16S_2018 <- learnErrors(filtRs_16S_2018, multithread = TRUE)

saveRDS(errF_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/errF_16S_2018")
saveRDS(errR_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/errR_16S_2018")
#The error rates for each possible nucleotide transition are shown. Points = error rates for each quality score, black line = estimated error rates from the machine-learning algorithm. Red lines shows expected error rate under nominal Q-score. If points fit well to the black line and error rate drops with increased quality score you can assume a good error rate.
Plot_Error_FWD_2018 <- plotErrors(errF_16S_2018, nominalQ = TRUE)
Plot_Error_REV_2018 <- plotErrors(errR_16S_2018, nominalQ = TRUE)
Plot_Error_FWD_2018
Plot_Error_REV_2018

saveRDS(Plot_Error_FWD_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/Plot_Error_FWD_2018")
saveRDS(Plot_Error_REV_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/Plot_Error_REV_2018")

#SECTION 4: Sample inference
#Apply the core sample inference algorithm to the filtered and trimmed sequence data
#Info: https://www.nature.com/articles/nmeth.3869#methods

dadaFs_16S_2018 <- dada(filtFs_16S_2018, err = errF_16S_2018, multithread = FALSE)
dadaRs_16S_2018 <- dada(filtRs_16S_2018, err = errR_16S_2018, multithread = FALSE)
dadaFs_16S_2018[[1]]
dadaRs_16S_2018[[1]]

saveRDS(dadaFs_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/dadaFs_16S_2018")
saveRDS(dadaRs_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/dadaRs_16S_2018")

mergers_16S_2018 <- mergePairs(dadaFs_16S_2018, filtFs_16S_2018, dadaRs_16S_2018, filtRs_16S_2018, verbose = TRUE)

# Inspect the merger data.frame from the first sample
head(mergers_16S_2018[[1]])

saveRDS(mergers_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/mergers_16S_2018")


#CONSTRUCT SEQUENCE TABLE
seqtab_16S_2018 <- makeSequenceTable(mergers_16S_2018)
dim(seqtab_16S_2018)
saveRDS(seqtab_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/seqtab_16S_2018")

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_16S_2018)))

# Remove all sequences outside the length range of 250-256.
seqtab_16S_2018_trim <- seqtab_16S_2018[,nchar(colnames(seqtab_16S_2018)) %in% seq(250,256)]
table(nchar(getSequences(seqtab_16S_2018_trim)))

saveRDS(seqtab_16S_2018_trim, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/seqtab_16S_2018_filtered")
seqtab_16S_2018_trim = readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/seqtab_16S_2018_filtered")

#REMOVE CHIMERAS
seqtab_nochim_16S_2018 <- removeBimeraDenovo(seqtab_16S_2018_trim, method = "consensus", multithread = FALSE, verbose = TRUE)

saveRDS(seqtab_nochim_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/seqtab_nochim_16S_2018")

# Inspect distribution of sequence lengths
seq_length_nochim_16S_2018 <- table(nchar(getSequences(seqtab_nochim_16S_2018)))
seq_length_nochim_16S_2018

write.csv(seq_length_nochim_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2016/seq_nochim_length_distribution")

dim(seqtab_16S_2018)
ncol(seqtab_16S_2018_trim)
ncol(seqtab_nochim_16S_2018)
dim(seqtab_nochim_16S_2018)
sum(seqtab_nochim_16S_2018) / sum(seqtab_16S_2018_trim)

#SECTION 5: Track reads through the pipeline
getN <- function(x)
  sum(getUniques(x))
getN

track_16S_2018 <-
  cbind(
    out_16S_2018,
    sapply(dadaFs_16S_2018, getN),
    sapply(dadaRs_16S_2018, getN),
    sapply(mergers_16S_2018, getN),
    rowSums(seqtab_nochim_16S_2018)
  )
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_16S_2018) <-
  c("input",
    "filtered",
    "denoisedF",
    "denoisedR",
    "merged",
    "nonchim")
rownames(track_16S_2018) <- sample.names_16S_2018
rownames(track_16S_2018)
print(track_16S_2018)
write.csv(track_16S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/16S_2018/track_16S_2018.csv")

#Outside of filtering, there should no step in which a majority of reads are lost.
#If a majority of reads failed to merge, you may need to revisit the truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon. #If a majority of reads were removed as chimeric, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification


#SECTION 6 - BRING IT ALL TOGETHER:
#IF YOU HAVE MULTIPLE SEQUENCE RUNS THIS IS THE TIME YOU MERGE THEM TOGETHER BEFORE ASSIGNING TAXONOMY, IF YOU'RE JUST WORKING WITH ONE SEQUENCE RUN IGNORE THIS STEP
row.names(seqtab_nochim_16S_2016)
row.names(seqtab_nochim_16S_2018) #Check you've got ALL your samples and no replicate names
st_16S <- mergeSequenceTables(seqtab_nochim_16S_2016, seqtab_nochim_16S_2018)
row.names(st_16S)
ncol(st_16S)

saveRDS(st_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/st_16S")

#SECTION 7 - Assign Taxonomy:
#You need to choose your reference database and have a file in your WD for this step
#SILVA and PR2 datasets available via DADA2 : https://benjjneb.github.io/dada2/training.html
merged_taxa_16S <- assignTaxonomy(st_16S, "tax/silva_nr99_v138.1_train_set.fa.gz", multithread = FALSE)
head(merged_taxa_16S)

saveRDS(merged_taxa_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/merged_16S_taxa")

#Optional: Exact species assignment
#This step adds exact species annotation to the fasta file
#Recent analysis suggest that exact 100% identity is the only appropriate way to assign species to 16S gene fragments, species-assignment training fastas are available for the SILVA 16S database
#Info: https://academic.oup.com/bioinformatics/article/34/14/2371/4913809?login=true

merged_taxa_16S <- addSpecies(merged_taxa_16S, "tax/silva_species_assignment_v138.1.fa.gz")
saveRDS(merged_taxa, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/merged_16S_taxa_species")

# I did this on a server as it was too large a file to process locally
merged_taxa_species <- readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/merged_taxa_16S_species")
merged_taxa_species

merged_taxa.print <-  merged_taxa_species #removing sequence rownames for display only
rownames(merged_taxa.print) <- NULL
head(merged_taxa.print)
write.table(merged_taxa.print, "merged_taxa_species_print_TEST.txt")
#Test your mock community
unqs.mock <- st_16S["Mock", ]
unqs.mock <-
  sort(unqs.mock[unqs.mock > 0], decreasing = TRUE) # Drop ASVs absent in the Mock
cat(
  "DADA2 inferred",
  length(unqs.mock),
  "sample sequences present in the Mock community.\n"
)
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <-
  sum(sapply(names(unqs.mock), function(x)
    any(grepl(x, mock.ref))))
cat("Of those,",
    sum(match.ref),
    "were exact matches to the expected reference sequences.\n") #this might say 0 due to the exact species assignment
#These can be pruned in phyloseq later

#9.1. LOAD METADATA FILE ####
metadata_16S <- read.table("canada_metadata_16S.txt", row.names = 1, header = TRUE)

#check that your samples and metadata table match (double check for typos to avoid much frustration down the line...)
rownames(metadata_16S) 
rownames(st_16S) #check your mock is still in there

#9.1. LOAD YOUR ASVs INTO PHYLOSEQ####
#Make a phyloseq object from your sequence table, tax file, and metadata. Phyloseq is a nice package for looking at amplicon data)
ps_16S <- phyloseq(
  otu_table(st_16S, taxa_are_rows = FALSE),
  sample_data(metadata_16S),
  tax_table(merged_taxa_species)
)

#Remove mock community from your phyloseq table
ps_16S <- prune_samples(sample_names(ps_16S) != "Mock", ps_16S)
ps_16S
sample_names(ps_16S)

#Move your DNA sequence data to a reference sequence slot in phyloseq for future reference
dna_16S <- Biostrings::DNAStringSet(taxa_names(ps_16S))
names(dna_16S) <- taxa_names(ps_16S)
ps_16S <- merge_phyloseq(ps_16S, dna_16S)
taxa_names(ps_16S) <- paste0("ASV", seq(ntaxa(ps_16S)))
ps_16S

#Access this data:
refseq(ps_16S)

#EXPORT RAW DATA####
#Export your raw ASV sequences as a FASTA file
ASV_seqs_16S <- refseq(ps_16S)
write.table(ASV_seqs_16S, file = "RAW_16S_ASVs.fasta", sep = "\t", quote=F, col.names=NA)

#Export the raw dada2 output as a merged ASV taxonomy table
#export an ASV count table
ps_asv_count_16S <- abundances(ps_16S)
ps_asv_count_16S <- as.data.frame(ps_asv_count_16S)
ps_asv_count_16S <- ps_asv_count_16S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_count_16S, file = "RAW_ps_ASV_count_16S.txt", sep = "\t", quote=F, col.names=NA)

#export a proportional abundance table of ASVs
ps_asv_ra_16S <- abundances(ps_16S, "compositional")
ps_asv_ra_16S <- as.data.frame(ps_asv_ra_16S)
ps_asv_ra_16S <- ps_asv_ra_16S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_ra_16S, file = "RAW_ps_ASV_proportional_count_16S.txt", sep = "\t", quote=F, col.names=NA)

#export a tax table
ps_tax <- tax_table(ps_16S)
ps_tax <- as.data.frame(ps_tax)
ps_tax <- ps_tax %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_tax, file = "RAW_ps_16S_taxa.txt", sep = "\t", quote=F, col.names=NA)

#export a merged ASV count and taxonomic assignment table
ASV_TAX_table_16S <- merge(ps_asv_count_16S, ps_tax, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_16S
write.table(ASV_TAX_table_16S, "RAW_merged_asv_tax_table_16S.txt", sep = "\t", quote=, col.names=NA)

#extract read count data and add to your metadata file
read_count <- readcount(ps_16S)
sample_data(ps_16S)$read_count <- read_count
sample_data(ps_16S)

#export your final metadata file as a .txt file
ps_meta <- meta(ps_16S)
row.names(ps_meta)
#Write it to a .txt file in your wd
write.table(ps_meta, file = "RAW_ps_metadata_16S.txt", sep = "\t", quote=F, col.names=NA)

#save raw (unfiltered) versions of your ASV files
saveRDS(dna_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/RAW_dna_16S")
saveRDS(ps_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/RAW_ps_16S")

ps_16S = readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/RAW_ps_16S")
summary_ps_16S <- as.data.frame(sort(sample_sums(ps_16S)))
summary_ps_16S 
summary(summary_ps_16S)

#Filtering ASV count data ####
#Investigate the total reads per sample and their distribution using a histogram
readsumsdf = data.frame(nreads = sort(taxa_sums(ps_16S), TRUE), sorted = 1:ntaxa(ps_16S), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps_16S), 
                                                        TRUE), sorted = 1:nsamples(ps_16S), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#Create table, number of features for each phyla
table(tax_table(ps_16S)[, "Phylum"], exclude = NULL)

# Abditibacteriota              Acidobacteriota             Actinobacteriota 
#42                         1007                         1286 
#Aenigmarchaeota                Altiarchaeota               Armatimonadota 
#3                            6                          462 
#Bacteroidota             Bdellovibrionota                Caldisericota 
#3195                         1483                           11 
#Campylobacterota                  Chloroflexi               Cloacimonadota 
#16                         2118                           17 
#Crenarchaeota                Cyanobacteria              Deferrisomatota 
#14                         2279                           11 
#Deinococcota                 Dependentiae             Desulfobacterota 
#36                          566                          444 
#Elusimicrobiota            Entotheonellaeota                Euryarchaeota 
#165                            1                           17 
#FCPU426            Fermentibacterota               Fibrobacterota 
#36                            1                          123 
#Firestonebacteria                   Firmicutes               Fusobacteriota 
#5                          656                            7 
#FW113              Gemmatimonadota                Halobacterota 
#2                          221                           47 
#Hydrogenedentes                Iainarchaeota             Latescibacterota 
#37                            4                           84 
#LCP-89             Margulisbacteria                       MBNT15 
#4                           24                           43 
#Methylomirabilota                Micrarchaeota                  Myxococcota 
#14                           11                         1425 
#Nanoarchaeota                        NB1-j                 Nitrospinota 
#308                           27                            7 
#Nitrospirota                        NKB15              Patescibacteria 
#67                            5                         2001 
#Planctomycetota               Proteobacteria                      RCP2-54 
#4675                         9947                            9 
#SAR324 clade(Marine group B)                Spirochaetota                  Sumerlaeota 
#525                          232                          169 
#Sva0485                         TA06             Thermoplasmatota 
#22                            1                           13 
#Verrucomicrobiota                        WOR-1                        WPS-2 
#3266                           14                           47 
#WS1                          WS2                          WS4 
#20                           14                           46 
#Zixibacteria                         <NA> 
#  26                         1267 #

#9.3. PREVALENCE BASED FILTERING####
#This is an unsupervised approach that relies on data in this experiment and a set parameter for filtering your ASVs, i.e. not based on taxonomic assignment

#Create a dataframe of ASV prevalence
#Prevalence defines the number of samples in which a taxon appears at once
prevdf_16S = apply(X = otu_table(ps_16S),
                       MARGIN = ifelse(taxa_are_rows(ps_16S), yes = 1, no = 2),
                       FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_16S = data.frame(Prevalence = prevdf_16S,
                            TotalAbundance = taxa_sums(ps_16S),
                            tax_table(ps_16S))
head(prevdf_16S)

write.csv(prevdf_16S, "prevalence_table_16S_raw.csv")
saveRDS(prevdf_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/prevdf_16S_raw")

#Are Divisions comprised of mostly low-prevalence features? compute the (1) average (mean) and (2) total (sum) prevalences of the features in each phylum
plyr::ddply(prevdf_16S, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

prevdf1 = subset(prevdf_16S, Phylum %in% get_taxa_unique(ps_16S, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps_16S),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

saveRDS(prevdf1, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/prevdf1_16S_raw")

prev_threshold <- 0.02 * nsamples(ps_16S) 
keepTaxa <- rownames(prevdf1)[prevdf1$Prevalence >= prev_threshold]
ps_16S_2_percent <- prune_taxa(keepTaxa, ps_16S)
ps_16S_2_percent

table(tax_table(ps_16S_2_percent)[,"Phylum"], exclude = NULL)

#Now look at how prevalence-based filtering has changed overall distribution of ASVs
prevdf <- apply(X = otu_table(ps_16S_2_percent), MARGIN = ifelse(taxa_are_rows(ps_16S_2_percent), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps_16S_2_percent), tax_table(ps_16S_2_percent))
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(ps_16S_2_percent, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps_16S_2_percent),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

saveRDS(ps_16S_2_percent, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/ps_16S_2_percent")
ps_16S_2_percent

#Make ps_16S_2_percent your main phyloseq file for ease
ps_16S = ps_16S_2_percent

#TAXONOMIC FILTERING####
# Create table, number of features for each phyla
table(tax_table(ps_16S)[, "Phylum"], exclude = NULL)
#Abditibacteriota              Acidobacteriota             Actinobacteriota 
#16                          567                          775 
#Aenigmarchaeota                Altiarchaeota               Armatimonadota 
#2                            2                          228 
#Bacteroidota             Bdellovibrionota                Caldisericota 
#1782                          510                            3 
#Campylobacterota                  Chloroflexi               Cloacimonadota 
#6                         1217                            6 
#Crenarchaeota                Cyanobacteria              Deferrisomatota 
#6                         1386                            5 
#Deinococcota                 Dependentiae             Desulfobacterota 
#14                          245                          253 
#Elusimicrobiota                Euryarchaeota                      FCPU426 
#45                            9                           15 
#Fermentibacterota               Fibrobacterota            Firestonebacteria 
#1                           39                            4 
#Firmicutes               Fusobacteriota                        FW113 
#254                            2                            1 
#Gemmatimonadota                Halobacterota              Hydrogenedentes 
#143                           30                           14 
#Iainarchaeota             Latescibacterota                       LCP-89 
#1                           43                            1 
#Margulisbacteria                       MBNT15            Methylomirabilota 
#8                           26                            6 
#Micrarchaeota                  Myxococcota                Nanoarchaeota 
#3                          741                           85 
#NB1-j                 Nitrospinota                 Nitrospirota 
#20                            3                           34 
#NKB15              Patescibacteria              Planctomycetota 
#2                          554                         2777 
#Proteobacteria                      RCP2-54 SAR324 clade(Marine group B) 
#4986                            2                          231 
#Spirochaetota                  Sumerlaeota                      Sva0485 
#112                           83                           11 
#TA06             Thermoplasmatota            Verrucomicrobiota 
#1                            5                         1677 
#WOR-1                        WPS-2                          WS1 
#2                           28                            8 
#WS2                          WS4                 Zixibacteria 
#4                           17                           10 
#<NA> 
#  420 

#Remove all eukaryotes
filterKingdoms = c("Eukaryota","NA") #change to your specific Supergroup(s)
ps_16S = subset_taxa(ps_16S, !Kingdom %in% filterKingdoms)
table(tax_table(ps_16S)[, "Kingdom"], exclude = NULL)
ps_16S
# Archaea Bacteria     <NA> 
#143    19328       10 

#Remove all NAs at the Phylum level
ps_16S <- subset_taxa(ps_16S, !is.na(Phylum) & !Phylum %in% c("", "NA"))
ps_16S

#Check the NAs have been removed and look at remaining groups, remove any with prevalence <3
table(tax_table(ps_16S)[, "Phylum"], exclude = NULL)

#Define phyla to filter such as anything that appears in a very low percent of all samples or plastid and mito
filterPhylum = c("Aenigmarchaeota", "Altiarchaeota", "Fermentibacterota","FW113","Fusobacteriota","Iainarchaeota","LCP-89","NKB15","RCP2-54","TA06","WOR-1") #change to your specific Supergroup(s)

# Remove filtered Phylum
ps_16S = subset_taxa(ps_16S, !Phylum %in% filterPhylum)
ps_16S

table(tax_table(ps_16S)[,"Phylum"], exclude = NULL) #check that <3 groups have been removed

#9.5. REMOVE UNWANTED SAMPLES ####
#ASV count per sample 
readcount(ps_16S)
#This shows that all samples have a good coverage and relatively similar number of reads, therefore none will be pruned

#9.6 FILTER OUT MITOCHONDRIA AND CHLOROPLAST
filter_chloro = c("Chloroplast") #change to your specific Supergroup(s)

# Remove chloro
ps_16S = subset_taxa(ps_16S, !Order %in% filter_chloro)
ps_16S

filter_mito = c("Mitochondria")

# Remove mito
ps_16S = subset_taxa(ps_16S, !Family %in% filter_mito)
ps_16S

#EXPORT FILTERED DATA####
#remove your DNA files again (you've already saved your raw data backup above)
dna_16S <- Biostrings::DNAStringSet(taxa_names(ps_16S))
names(dna_16S) <- taxa_names(ps_16S)
ps_16S <- merge_phyloseq(ps_16S, dna_16S)
taxa_names(ps_16S) <- paste0("ASV", seq(ntaxa(ps_16S)))
ps_16S
#Access this data:
refseq(ps_16S)

#9.5. SAVE YOUR FILTERED PHYLOSEQ OBJECT
saveRDS(dna_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/filt_dna_16S")
saveRDS(ps_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/filt_ps_16S")

#Export your FILTERED ASV sequences as a FASTA file
ASV_seqs_16S <- refseq(ps_16S)
write.table(ASV_seqs_16S, file = "FILT_16S_ASVs.fasta", sep = "\t", quote=F, col.names=NA)

#Export the raw dada2 output as a merged ASV taxonomy table
#export an ASV count table
ps_asv_count_16S <- abundances(ps_16S)
ps_asv_count_16S <- as.data.frame(ps_asv_count_16S)
ps_asv_count_16S <- ps_asv_count_16S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_count_16S, file = "FILT_ps_ASV_count_16S.txt", sep = "\t", quote=F, col.names=NA)

#export a proportional abundance table of ASVs
ps_asv_ra_16S <- abundances(ps_16S, "compositional")
ps_asv_ra_16S <- as.data.frame(ps_asv_ra_16S)
ps_asv_ra_16S <- ps_asv_ra_16S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_ra_16S, file = "FILT_ps_ASV_proportional_count_16S.txt", sep = "\t", quote=F, col.names=NA)

#export a tax table
ps_tax <- tax_table(ps_16S)
ps_tax <- as.data.frame(ps_tax)
ps_tax <- ps_tax %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_tax, file = "FILT_ps_16S_taxa.txt", sep = "\t", quote=F, col.names=NA)

#export a merged ASV count and taxonomic assignment table
ASV_TAX_table_16S <- merge(ps_asv_count_16S, ps_tax, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_16S
write.table(ASV_TAX_table_16S, "FILT_merged_asv_tax_table_16S.txt", sep = "\t", quote=, col.names=NA)

#export a merged RA ASV count and taxonomic assingment table
ASV_TAX_table_16S_ra <- merge(ps_asv_ra_16S, ps_tax, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_16S_ra
write.table(ASV_TAX_table_16S_ra, "FILT_merged_asv_tax_table_16S_RA.txt", sep = "\t", quote=, col.names=NA)

#extract read count data and add to your metadata file
read_count <- readcount(ps_16S)
sample_data(ps_16S)$read_count <- read_count
sample_data(ps_16S)

#export your final metadata file as a .txt file
ps_meta <- meta(ps_16S)
row.names(ps_meta)
#Write it to a .txt file in your wd
write.table(ps_meta, file = "FILT_ps_metadata_16S.txt", sep = "\t", quote=F, col.names=NA)

#save raw (unfiltered) versions of your ASV files
saveRDS(dna_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FILT_dna_16S")
saveRDS(ps_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FILT_ps_16S")

ps_16S

#ADD PHYLOGENETIC FILE
ps.merged = readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/ps_merged_16S_2_4_22")
ps.merged

#output:
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 17952 taxa and 97 samples ]
#sample_data() Sample Data:       [ 97 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 17952 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 17952 tips and 17950 internal nodes ]
#refseq()      DNAStringSet:      [ 17952 reference sequences ]

#9.9. PHYLOGENETIC AGGLOMERATION####
#This step uses the tip_glom feature of phyloseq to cluster ASVs based on their phylogenetic relatedness
#This also had to be carried out on a server using the script "glom_script.R"

#I tried this at 0.2, 0.4, 0.02, and 0.01 and selected 0.01 to match 18S
h1 = 0.01

ps_glom_0.01 <- readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/ps_merged_16S_glom_0.01_4_4")
ps_glom_0.01

#output
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 15250 taxa and 97 samples ]
#sample_data() Sample Data:       [ 97 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 15250 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 15250 tips and 15249 internal nodes ]
#refseq()      DNAStringSet:      [ 15250 reference sequences ]
#plot_tree(ps_glom_0.02, color="region", label.tips="Genus")

#Look at how agglomeration changes the phylogenetic tree
multiPlotTitleTextSize = 15
test_counts = plot_tree(ps.merged,
                        method = "treeonly",
                        ladderize = "left",
                        title = "16S Before Agglomeration") + theme(plot.title = element_text(size = multiPlotTitleTextSize))
tipglom_counts_0.01 = plot_tree(ps_glom_0.01,
                               method = "treeonly",
                               ladderize = "left",
                               title = "16S_Glom_0.01") + theme(plot.title = element_text(size = multiPlotTitleTextSize))
grid.arrange(nrow = 1, test_counts, tipglom_counts_0.01)

#START HERE IN MORNING - CHANGE TO 16s
#9.10. EXPORT YOUR FINAL PHYLOGENETIC AGGLOMERATED DATA TABLE
#export an ASV count table
ps_asv_count_16S <- abundances(ps_glom_0.01)
ps_asv_count_16S <- as.data.frame(ps_asv_count_16S)
ps_asv_count_16S <- ps_asv_count_16S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_count_16S, file = "FINAL_ps_ASV_count_16S_glom_0_01.txt", sep = "\t", quote=F, col.names=NA)

#export a proportional abundance table of ASVs
ps_asv_ra_16S <- abundances(ps_glom_0.01,"compositional")
ps_asv_ra_16S <- as.data.frame(ps_asv_ra_16S)
ps_asv_ra_16S <- ps_asv_ra_16S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_ra_16S, file = "FINAL_ps_ASV_proportional_count_16S_glom_0_01.txt", sep = "\t", quote=F, col.names=NA)

#export a tax table
ps_tax_16S <- tax_table(ps_glom_0.01)
ps_tax_16S <- as.data.frame(ps_tax_16S)
ps_tax_16S <- ps_tax_16S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_tax_16S, file = "FINAL_ps_16S_taxa_glom_0_01.txt", sep = "\t", quote=F, col.names=NA)

#export a merged ASV count and taxonomic assignment table
ASV_TAX_table_16S <- merge(ps_asv_count_16S, ps_tax_16S, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_16S
write.table(ASV_TAX_table_16S, "FINAL_merged_asv_tax_table_16S_glom_0_01.txt", sep = "\t", quote=, col.names=NA)

#export a merged proportional (relative abundance) ASV count and taxonomic assignment table
ASV_TAX_table_16S_RA <- merge(ps_asv_ra_16S, ps_tax_16S, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_16S_RA
write.table(ASV_TAX_table_16S_RA, "FINAL_merged_asv_tax_table_16S_glom_0_01_RA.txt", sep = "\t", quote=, col.names=NA)

#export your ASV sequences as a FASTA file
ASV_seqs_16S <- refseq(ps_glom_0.01)
write.table(ASV_seqs_16S, file = "FINAL_16S_ASVs_glom_0_01.fasta", sep = "\t", quote=F, col.names=NA)

#export read count
read_count <- readcount(ps_glom_0.01)

#Add this to your final phyloseq metadata object
sample_data(ps_glom_0.01)$read_count <- read_count
sample_data(ps_glom_0.01)

#psmelt function allows you to convert phyloseq objects into R dataframes, useful for future work with ggplot
meta <- sample_data(ps_glom_0.01)
write.table(meta, "FINAL_ps_metadata_16S_glom_0_01.txt", sep = "\t", quote=, col.names=NA)

#Save your final version of your phyloseq object for easy reloading later on
ps_glom_0.01
saveRDS(ps_glom_0.01, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FINAL_ps_16S")

#Create a histogram of your final ASV and Sample counts
readsumsdf = data.frame(nreads = sort(taxa_sums(ps_glom_0.01), TRUE), sorted = 1:ntaxa(ps_glom_0.01), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps_glom_0.01), 
                                                        TRUE), sorted = 1:nsamples(ps_glom_0.01), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#Create a plot of the final prevalence of ASVs faceted taxonomically by Division 
prevdf_16S = apply(X = otu_table(ps_glom_0.01),
                   MARGIN = ifelse(taxa_are_rows(ps_glom_0.01), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_16S = data.frame(Prevalence = prevdf_16S,
                        TotalAbundance = taxa_sums(ps_glom_0.01),
                        tax_table(ps_glom_0.01))
head(prevdf_16S)

write.table(prevdf_16S, "FINAL_prevalence_table_16S.txt")
saveRDS(prevdf_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FINAL_prevdf_16S")

#Are Divisions comprised of mostly low-prevalence features? compute the (1) average (mean) and (2) total (sum) prevalences of the features in each phylum
plyr::ddply(prevdf_16S, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
prevdf1 <- subset(prevdf_16S, Phylum %in% get_taxa_unique(ps_glom_0.01, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps_glom_0.01),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#END OF SCRIPT####
