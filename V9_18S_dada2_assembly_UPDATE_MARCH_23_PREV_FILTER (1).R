#Assembly of V9 18S sequence data using dada2 ASV assembly 
##Author: Patrick M. Hooper
##Github:https://github.com/paddyhooper/dada2
##Created: 01.04.22

#Along the way I used saveRDS to save R products, this can save a lot of time in the future as you can just ##Reload previously run command outputs using the script in my GitHub page: "Reloading your DADA2 analyses"

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

##YOU'RE READY TO BEGIN WITH DADA2
#My sequencing files are from two separate sequencing runs, DADA2 requires you to assemble sequencing runs separately and merge them after assembly
#I have two sequencing runs from 2016 and 2018, I started with the 2016 sequences

#1.1. IDENTIFY PATH TO YOUR 18S FASTA SEQUENCING FILES####
path <- "C:/Users/pmh36/Documents/2016_18S_Files"
list.files(path) # check your files are all present

#Forward and reverse fastq file names have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs_18S_2016 <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs_18S_2016 <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
fnFs_18S_2016 # fwd reads
fnRs_18S_2016 # rev reads

#1.2. EXTRACT SAMPLE NAMES FROM PATH FILES####
#Assumes file names have format: SAMPLENAME_XXX.fastq
sample.names_18S_2016 <- sapply(strsplit(basename(fnFs_18S_2016), "_"), `[`, 1)
sample.names_18S_2016 #double check you have all your samples

#1.3. INSPECT READ QUALITY####
#Grey scale heat map is the frequency of each quality score at each base position
#Mean quality score = green line
#Quartiles of quality score = orange line
#Red line = scaled proportion of reads extend to at least that position, not useful for Illumina as all reads same length
FWD_Quality_Profile_2016 <- plotQualityProfile(fnFs_18S_2016[1:4])
REV_Quality_Profile_2016 <- plotQualityProfile(fnRs_18S_2016[1:4])
FWD_Quality_Profile_2016
REV_Quality_Profile_2016

saveRDS(FWD_Quality_Profile_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/FWD_Quality_Profile_2016")
saveRDS(REV_Quality_Profile_2016 , file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/REV_Quality_Profile_2016")

#2.1. FILTER AND TRIM####
# Place filtered files in 'filtered/' subdirectory
filtFs_18S_2016 <- file.path(path, "filtered", paste0(sample.names_18S_2016, "_F_filt.fastq.gz"))
filtRs_18S_2016 <- file.path(path, "filtered", paste0(sample.names_18S_2016, "_R_filt.fastq.gz"))
names(filtFs_18S_2016) <- sample.names_18S_2016
names(filtRs_18S_2016) <- sample.names_18S_2016

saveRDS(filtFs_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/filtFs_18S_2016")
saveRDS(filtRs_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/filtRs_18S_2016")

out_18S_2016 <- filterAndTrim(
  fnFs_18S_2016,
  filtFs_18S_2016,
  fnRs_18S_2016,
  filtRs_18S_2016,
  trimLeft = 5,
  truncLen = c(100, 100),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE
) # On Windows set multithread = FALSE

saveRDS(out_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/out_18S_2016")

#3.1. LEARN ERROR RATES####
errF_18S_2016 <- learnErrors(filtFs_18S_2016, multithread = TRUE)
errR_18S_2016 <- learnErrors(filtRs_18S_2016, multithread = TRUE)
errF_18S_2016
errR_18S_2016

saveRDS(errF_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/errF_18S_2016")
saveRDS(errR_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/errR_18S_2016")

#3.2. Plot your error graphs
#The error rates for each possible nucleotide transition are shown. Points = error rates for each quality score, black line = estimated error rates from the machine-learning algorithm. Red lines shows expected error rate under nominal Q-score.
#INTERPRETATION: If points fit well to the black line and error rate drops with increased quality score you can assume a good error rate.

Plot_Error_FWD <- plotErrors(errF_18S_2016, nominalQ=TRUE)#FWD read errors
Plot_Error_REV <- plotErrors(errR_18S_2016, nominalQ=TRUE) #REV read errors
Plot_Error_FWD
Plot_Error_REV

saveRDS(Plot_Error_FWD, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/Plot_Error_FWD")
saveRDS(Plot_Error_REV, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/Plot_Error_REV")

#4.1. SAMPLE INFERENCE####
#Apply the core sample inference algorithm to the filtered and trimmed sequence data
#Info: https://www.nature.com/articles/nmeth.3869#methods
dadaFs_18S_2016 <- dada(filtFs_18S_2016, err = errF_18S_2016, multithread = FALSE)
dadaRs_18S_2016 <- dada(filtRs_18S_2016, err = errR_18S_2016, multithread = FALSE)

saveRDS(dadaFs_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/dadaFs_18S_2016")
saveRDS(dadaRs_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/dadaRs_18S_2016")

#Inspect the first line
dadaFs_18S_2016[[1]]

#4.2. Merge your filtered sequence files with your sample inference
mergers_18S_2016 <-
  mergePairs(dadaFs_18S_2016,
             filtFs_18S_2016,
             dadaRs_18S_2016,
             filtRs_18S_2016,
             verbose = TRUE)

#4.3. Inspect the merger data.frame from the first sample
head(mergers_18S_2016[[1]])

saveRDS(mergers_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/mergers_18S_2016")

#5.1. CONSTRUCT SEQUENCE TABLE####
seqtab_18S_2016 <- makeSequenceTable(mergers_18S_2016)
dim(seqtab_18S_2016)

saveRDS(seqtab_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/seqtab_18S_2016")

#5.2. Inspect distribution of sequence lengths
seq_length_18S_2016 = table(nchar(getSequences(seqtab_18S_2016)))
seq_length_18S_2016
write.csv(seq_length_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/seq_length_distribution")

#5.3 Trim those sequences that are significantly larger than the expected primer length ~130 for V9 
seqtab_18S_2016_trim<- seqtab_18S_2016[,nchar(colnames(seqtab_18S_2016)) %in% seq(95,150)]
table(nchar(getSequences(seqtab_18S_2016_trim)))

saveRDS(seqtab_18S_2016_trim, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/seqtab_18S_2016_new_trim_filtered")

#5.3. REMOVE CHIMERAS
seqtab_nochim_18S_2016 <- removeBimeraDenovo(seqtab_18S_2016_trim, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab_nochim_18S_2016)
table(nchar(getSequences(seqtab_nochim_18S_2016)))

saveRDS(seqtab_nochim_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/seqtab_nochim_18S_2016")

#How many sequences were retained after chimeras removed
ncol(seqtab_18S_2016)
ncol(seqtab_nochim_18S_2016)
sum(seqtab_nochim_18S_2016) / sum(seqtab_18S_2016) *100 # % of non-chimeras

#6.1. TRACK YOUR READS THROUGH THE PIPELINE####
getN <- function(x)
  sum(getUniques(x))
getN

track_18S_2016 <-
  cbind(
    out_18S_2016,
    sapply(dadaFs_18S_2016, getN),
    sapply(dadaRs_18S_2016, getN),
    sapply(mergers_18S_2016, getN),
    rowSums(seqtab_nochim_18S_2016)
  )

#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track_18S_2016) <-
  c("input",
    "filtered",
    "denoisedF",
    "denoisedR",
    "merged",
    "nonchim")
rownames(track_18S_2016) <- sample.names_18S_2016
print(track_18S_2016)

saveRDS(track_18S_2016, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/track_18S_2016")

#Outside of filtering, there should no step in which a majority of reads are lost
#If a majority of reads failed to merge, you may need to revisit the truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon
#If a majority of reads were removed as chimeric, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification

#OPTIONAL 2: IDENTIFY PATH TO YOUR SECOND FOLDER OF SEQUENCING FILES####
#Skip if you're only analysing one sequencing run
path <- "C:/Users/pmh36/Documents/2018_18S_Files"
list.files(path) # check your files are all present

#Forward and reverse fastq file names have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs_18S_2018 <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs_18S_2018 <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
fnFs_18S_2018 # fwd reads
fnRs_18S_2018 # rev reads

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names_18S_2018 <- sapply(strsplit(basename(fnFs_18S_2018), "_"), `[`, 1)
sample.names_18S_2018
###If you need to remove primers in this dataset, repeat cut adapt steps above###

#SECTION 1: Inspect read quality profiles for FWD reads for first 4 sequences
FWD_Quality_Profile_2018 <- plotQualityProfile(fnFs_18S_2018[1:4])
REV_Quality_Profile_2018 <- plotQualityProfile(fnRs_18S_2018[1:4])
FWD_Quality_Profile_2018
REV_Quality_Profile_2018

saveRDS(FWD_Quality_Profile_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/FWD_Quality_Profile_2018_new_trim")
saveRDS(REV_Quality_Profile_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/REV_Quality_Profile_2018_new_trim")

#SECTION 2: Filter and Trim
# Place filtered files in filtered/ subdirectory
filtFs_18S_2018 <- file.path(path, "filtered", paste0(sample.names_18S_2018, "_F_filt.fastq.gz"))
filtRs_18S_2018 <- file.path(path, "filtered", paste0(sample.names_18S_2018, "_R_filt.fastq.gz"))
names(filtFs_18S_2018) <- sample.names_18S_2018
names(filtRs_18S_2018) <- sample.names_18S_2018

out_18S_2018 <- filterAndTrim(
  fnFs_18S_2018,
  filtFs_18S_2018,
  fnRs_18S_2018,
  filtRs_18S_2018,
  trimLeft = 5,
  truncLen = c(100, 100),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE
) # On Windows set multithread = FALSE

#save an RDS file
saveRDS(out_18S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/out_18S_2018_new_trim")

#SECTION 3: Learn the error rates
errF_18S_2018 <- learnErrors(filtFs_18S_2018, multithread = TRUE)
errR_18S_2018 <- learnErrors(filtRs_18S_2018, multithread = TRUE)
errF_18S_2018
errR_18S_2018

saveRDS(errF_18S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/errF_18S_2018_new_trim")
saveRDS(errR_18S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/errR_18S_2018_new_trim")

errF_18S_2018 <- readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/errF_18S_2018_new_trim")
errR_18S_2018 <- readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/errR_18S_2018_new_trim")

#The error rates for each possible nucleotide transition are shown. Points = error rates for each quality score, black line = estimated error rates from the machine-learning algorithm. Red lines shows expected error rate under nominal Q-score. If points fit well to the black line and error rate drops with increased quality score you can assume a good error rate.

Plot_Error_FWD_2018 <- plotErrors(errF_18S_2018, nominalQ=TRUE) #FWD read errors
Plot_Error_REV_2018 <- plotErrors(errR_18S_2018, nominalQ=TRUE) #REV read errors
Plot_Error_FWD_2018
Plot_Error_REV_2018

saveRDS(Plot_Error_FWD_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/Plot_Error_FWD_2018_new_trim")
saveRDS(Plot_Error_REV_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/Plot_Error_REV_2018_new_trim")

#SECTION 4: Sample inference
#Apply the core sample inference algorithm to the filtered and trimmed sequence data
#Info: https://www.nature.com/articles/nmeth.3869#methods
dadaFs_18S_2018 <- dada(filtFs_18S_2018, err=errF_18S_2018, multithread=FALSE)
saveRDS(dadaFs_18S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/dadaFs_18S_2018_new_trim")
dadaRs_18S_2018 <- dada(filtRs_18S_2018, err=errR_18S_2018, multithread=FALSE)
saveRDS(dadaRs_18S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/dadaRs_18S_2018_new_trim")
dadaFs_18S_2018[[1]]
dadaRs_18S_2018[[1]]

#SECTION 5: MERGE YOUR SAMPLE INFERENCE WITH YOUR SEQUENCE READS
mergers_18S_2018 <- mergePairs(dadaFs_18S_2018, filtFs_18S_2018, dadaRs_18S_2018, filtRs_18S_2018, verbose=TRUE)
head(mergers_18S_2018)

saveRDS(mergers_18S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/mergers_18S_2018_new_trim")

#Inspect the merger data.frame from the first sample
head(mergers_18S_2018[[1]])

#SECTION 6: CONSTRUCT SEQUENCE TABLE
seqtab_18S_2018 <- makeSequenceTable(mergers_18S_2018)
dim(seqtab_18S_2018)

saveRDS(seqtab_18S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/seqtab_18S_2018_new_trim")

#5.2. Inspect distribution of sequence lengths
seq_length_18S_2018 = table(nchar(getSequences(seqtab_18S_2018)))
seq_length_18S_2018
write.csv(seq_length_18S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/seq_length_distribution")

#5.3 Trim those sequences that are significantly larger than the expected primer length ~130 for V9 
seqtab_18S_2018_trim<- seqtab_18S_2018[,nchar(colnames(seqtab_18S_2018)) %in% seq(95,150)]
table(nchar(getSequences(seqtab_18S_2018_trim)))

saveRDS(seqtab_18S_2018_trim, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/seqtab_18S_2018_new_trim_filtered")

#check sample names
rownames(seqtab_18S_2018)

#REMOVE CHIMERAS
seqtab_nochim_18S_2018 <- removeBimeraDenovo(seqtab_18S_2018_trim, method="consensus", multithread=FALSE, verbose=TRUE)

saveRDS(seqtab_nochim_18S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/seqtab_nochim_18S_2018_new_trim")

#change sample names
rownames(seqtab_nochim_18S_2018)
rownames(seqtab_nochim_18S_2018) <- sample.names_18S_2018
rownames(seqtab_nochim_18S_2018)
table(nchar(getSequences(seqtab_nochim_18S_2018)))

#Percentage of chimera sequences
ncol(seqtab_18S_2018) #how many ASVs 
ncol(seqtab_nochim_18S_2018) # how many ASVs
dim(seqtab_nochim_18S_2018)
rownames(seqtab_nochim_18S_2018)
sum(seqtab_nochim_18S_2018)/sum(seqtab_18S_2018)*100 # % of non-chimeras

#SECTION 7: Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
getN

track_18S_2018 <-
  cbind(
    out_18S_2018,
    sapply(dadaFs_18S_2018, getN),
    sapply(dadaRs_18S_2018, getN),
    sapply(mergers_18S_2018, getN),
    rowSums(seqtab_nochim_18S_2018)
  )

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_18S_2018) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_18S_2018) <- sample.names_18S_2018
rownames(track_18S_2018)
print(track_18S_2018)

saveRDS(track_18S_2018, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2018/track_18S_2018")
write.table(track_18S_2018, file = "track_18S_2018_new_trim.txt")

#IF YOU WERE WORKING ON ONE SEQUENCING RUN START AGAIN HERE#

#7.1. COMBINE YOUR SEQUENCING RUNS ####
#Check you've got ALL your samples and no replicate names
seqtab_nochim_18S_2016 = readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_2016/seqtab_nochim_18S_2016")
row.names(seqtab_nochim_18S_2016)
row.names(seqtab_nochim_18S_2018) 

#IF YOU HAVE MULTIPLE SEQUENCE RUNS THIS IS THE TIME YOU MERGE THEM TOGETHER BEFORE ASSIGNING TAXONOMY, IF YOU'RE JUST WORKING WITH ONE SEQUENCE RUN IGNORE THIS STEP
st_18S <- mergeSequenceTables(seqtab_nochim_18S_2016, seqtab_nochim_18S_2018)
head(st_18S)
row.names(st_18S)

# Inspect distribution of sequence lengths
table(nchar(getSequences(st_18S)))

saveRDS(st_18S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/st_18S_new_trim")#This is a big one, make sure you've saved this and you won't have to re-analyse the whole thing! #note this changes the file directory from the previous files

#8.1. ASSIGN TAXONOMY WITH PR2 ####
#PR2 didn't work locally so this has to be run on a server using the 'PR2 18S taxonomic classification.R'  script file and the RScript function on a UNIX server
#merged_taxa_18S <- assignTaxonomy(st_18S, "pr2_version_4.14.0_SSU_dada2.fasta.gz", multithread=TRUE, taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"))

merged_taxa_18S <- readRDS("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/merged_taxa_18S_31")
merged_taxa.print <- merged_taxa_18S #removing sequence rownames for display only
rownames(merged_taxa.print) <- NULL
head(merged_taxa.print)

#9.1. LOAD METADATA FILE ####
metadata_18S <- read.table("canada_metadata_18S.txt",	row.names=1, header = TRUE)
metadata_18S

#check that your samples and metadata table match (double check for typos to avoid much frustration down the line...)
rownames(metadata_18S)
rownames(st_18S) 

#9.1. LOAD YOUR ASVs INTO PHYLOSEQ####
#Make a phyloseq object from your sequence table, tax file, and metadata. Phyloseq is a nice package for looking at amplicon data)
ps_18S_new_trim <- phyloseq(
  otu_table(st_18S, taxa_are_rows = FALSE),
  #taxa_are_rows is FALSE if importing directly from dada2 but TRUE if importing from QIIME2
  sample_data(metadata_18S),
  tax_table(merged_taxa_18S)
)
ps_18S_new_trim #check your phyloseq object, how many ASVs and samples
sample_names(ps_18S_new_trim) 

#Move your DNA sequence data to a reference sequence slot in phyloseq for future reference
dna_18S_new_trim <- Biostrings::DNAStringSet(taxa_names(ps_18S_new_trim))
names(dna_18S_new_trim) <- taxa_names(ps_18S_new_trim)
ps_18S_new_trim <- merge_phyloseq(ps_18S_new_trim, dna_18S_new_trim)
taxa_names(ps_18S_new_trim) <- paste0("ASV", seq(ntaxa(ps_18S_new_trim)))
ps_18S_new_trim

#Access this data:
refseq(ps_18S_new_trim)

#EXPORT RAW DATA####
#Export your raw ASV sequences as a FASTA file
ASV_seqs_new_trim <- refseq(ps_18S_new_trim)
write.table(ASV_seqs_new_trim, file = "RAW_18S_ASVs_new_trim.fasta", sep = "\t", quote=F, col.names=NA)

#Export the raw dada2 output as a merged ASV taxonomy table
#export an ASV count table
ps_asv_count_18S_new_trim <- abundances(ps_18S_new_trim)
ps_asv_count_18S_new_trim <- as.data.frame(ps_asv_count_18S_new_trim)
ps_asv_count_18S_new_trim <- ps_asv_count_18S_new_trim %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_count_18S_new_trim, file = "RAW_ps_ASV_count_18S_new_trim.txt", sep = "\t", quote=F, col.names=NA)

#export a proportional abundance table of ASVs
ps_asv_ra_18S <- abundances(ps_18S_new_trim, "compositional")
ps_asv_ra_18S <- as.data.frame(ps_asv_ra_18S)
ps_asv_ra_18S <- ps_asv_ra_18S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_ra_18S, file = "RAW_ps_ASV_proportional_count_18S_new_trim.txt", sep = "\t", quote=F, col.names=NA)

#export a tax table
ps_tax <- tax_table(ps_18S_new_trim)
ps_tax <- as.data.frame(ps_tax)
ps_tax <- ps_tax %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_tax, file = "RAW_ps_18S_taxa_new_trim.txt", sep = "\t", quote=F, col.names=NA)

#export a merged ASV count and taxonomic assignment table
ASV_TAX_table_18S_new_trim <- merge(ps_asv_count_18S_new_trim, ps_tax, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_18S_new_trim
write.table(ASV_TAX_table_18S_new_trim, "RAW_merged_asv_tax_table_18S_new_trim.txt", sep = "\t", quote=, col.names=NA)

#extract read count data and add to your metadata file
read_count <- readcount(ps_18S_new_trim)
sample_data(ps_18S_new_trim)$read_count <- read_count
sample_data(ps_18S_new_trim)

#export your final metadata file as a .txt file
ps_meta_new_trim <- meta(ps_18S_new_trim)
row.names(ps_meta_new_trim)
#Write it to a .txt file in your wd
write.table(ps_meta_new_trim, file = "RAW_ps_metadata_new_trim.txt", sep = "\t", quote=F, col.names=NA)

#save raw (unfiltered) versions of your ASV files
saveRDS(dna_18S_new_trim, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/RAW_dna_18S_new_trim")
saveRDS(ps_18S_new_trim, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/RAW_ps_18S_new_trim")

#checking the prevalence trimmed version
ps_18S_new_trim <- readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/RAW_ps_18S_new_trim")

#check sequence depth
ps_18S = readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/RAW_ps_18S_new_trim")
summary_ps_18S <- as.data.frame(sort(sample_sums(ps_18S)))
summary_ps_18S 
summary(summary_ps_18S)

#9.2. FILTERING ASV TABLES IN PHYLOSEQ ####
#Phyloseq provides useful tools for filtering, subsetting, and agglomerating taxa
#We will now use phyloseq to edit the content of the 18S ASV sequence table before extracting the final standard products (count table, tax table, metadata)

#Investigate the total reads per sample and their distribution using a histogram
readsumsdf = data.frame(nreads = sort(taxa_sums(ps_18S_new_trim), TRUE), sorted = 1:ntaxa(ps_18S_new_trim), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps_18S_new_trim), 
                                                        TRUE), sorted = 1:nsamples(ps_18S_new_trim), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")

p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")


## i didn't do the below information

#9.3. PREVALENCE BASED FILTERING####
#This is an unsupervised approach that relies on data in this experiment and a set parameter for filtering your ASVs, i.e. not based on taxonomic assignment

#Create a dataframe of ASV prevalence
#Prevalence defines the number of samples in which a taxon appears at once
prevdf_18S_raw = apply(X = otu_table(ps_18S_new_trim),
                   MARGIN = ifelse(taxa_are_rows(ps_18S_new_trim), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_18S_raw = data.frame(Prevalence = prevdf_18S_raw,
                             TotalAbundance = taxa_sums(ps_18S_new_trim),
                             tax_table(ps_18S_new_trim))
head(prevdf_18S_raw)

write.csv(prevdf_18S_raw, "prevalence_table_18S_raw.csv")
saveRDS(prevdf_18S_raw, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/prevdf_18S_raw")

#Are Divisions comprised of mostly low-prevalence features? compute the (1) average (mean) and (2) total (sum) prevalences of the features in each phylum
plyr::ddply(prevdf_18S_raw, "Division", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

prevdf1 <- subset(prevdf_18S_raw, Division %in% get_taxa_unique(ps_18S_new_trim, "Division"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps_18S_new_trim),color=Division)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Division) + theme(legend.position="none")

prev_threshold <- 0.02 * nsamples(ps_18S_new_trim) 
keepTaxa <- rownames(prevdf1)[prevdf1$Prevalence >= prev_threshold]
ps_18S_raw_2_percent <- prune_taxa(keepTaxa, ps_18S_new_trim)
ps_18S_raw_2_percent
ps_18S_new_trim
table(tax_table(ps_18S_raw_2_percent)[,"Division"], exclude = NULL)
table(tax_table(ps_18S_new_trim)[,"Division"], exclude = NULL)

#9.4. TAXONOMIC FILTERING ####
rank_names(ps_18S_new_trim) # show available rank names in the dataset

# Create table, number of features for each taxonomic division
table(tax_table(ps_18S_raw_2_percent)[, "Division"], exclude = NULL)

#OUTPUT:
#Acidobacteria     Actinobacteria        Apicomplexa     Apusomonadidae          Breviatea 
#61                  1                163                  6                 10 
#Centroheliozoa           Cercozoa        Chloroflexi        Chlorophyta   Choanoflagellida 
#64                899                  1                733                 33 
#Chrompodellids         Ciliophora             Conosa        Cryptophyta      Cyanobacteria 
#10                857                142                 27                  2 
#Dependentiae     Dinoflagellata            Discoba Epsilonbacteraeota       Eukaryota_XX 
#27                363                440                  1                  8 
#Euryarchaeota         Firmicutes       Foraminifera              Fungi   Gemmatimonadetes 
#20                 30                 47               1974                 31 
#Glaucophyta         Haptophyta   Hemimastigophora        Hilomonadea    Hydrogenedentes 
#4                 22                  2                 17                  1 
#Kiritimatiellaeota             Lobosa      Mesomycetozoa         Metamonada            Metazoa 
#10                416                 54                 14               1366 
#Nitrospirae         Ochrophyta           Opalozoa    Patescibacteria            Picozoa 
#7               1054                 87                  5                  1 
#Planctomycetes     Proteobacteria        Pseudofungi         Rhodophyta          Sagenista 
#334               1482                306                 46                239 
#Stramenopiles_X       Streptophyta          Telonemia     Thaumarchaeota    Verrucomicrobia 
#1                 45                  2                  2                252 
#<NA> 
#7270 #quite a large amount of NA reads

#remove any bacteria ASVs
ps_18S_raw_2_percent = subset_taxa(ps_18S_raw_2_percent, Kingdom == "Eukaryota")
ps_18S_raw_2_percent

#remove anything that is considered ambiguous annotation, i.e. NA annotation at Supergroup level
ps_18S_raw_2_percent <- subset_taxa(ps_18S_raw_2_percent, !is.na(Supergroup) & !Supergroup %in% c("Eukaryota_X", "NA"))
ps_18S_raw_2_percent

#create a table of the number of ASVs in each Supergroup (L2)
table(tax_table(ps_18S_raw_2_percent)[,"Division"], exclude = NULL) 
#This shows us which divisions only feature a few ASVS, these might be worth filtering
#Apicomplexa   Apusomonadidae 
#163                6 
#Breviatea   Centroheliozoa 
#10               64 
#Cercozoa      Chlorophyta 
#899              733 
#Choanoflagellida   Chrompodellids 
#33               10 
#Ciliophora           Conosa 
#857              142 
#Cryptophyta   Dinoflagellata 
#27              363 
#Discoba     Foraminifera 
#440               47 
#Fungi      Glaucophyta 
#1974                4 
#Haptophyta      Hilomonadea 
#22               17 
#Lobosa    Mesomycetozoa 
#416               54 
#Metamonada          Metazoa 
#14             1366 
#Ochrophyta         Opalozoa 
#1054               87 
#Picozoa      Pseudofungi 
#1              306 
#Rhodophyta        Sagenista 
#46              239 
#Stramenopiles_X     Streptophyta 
#1               45 
#Telonemia       <NA> 
#  2             1026 

#Compute the prevalence of features in each phylum, how often do they occur?
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#Define phyla to filter such as anything that appears in a very low percent of all samples or plastid and mito
filterDivision = c("Stramenopiles_X", "Picozoa", "Telonemia") #change to your specific Supergroup(s)

# Remove filtered Supergroups
ps_18S_raw_2_percent = subset_taxa(ps_18S_raw_2_percent, !Division %in% filterDivision)

table(tax_table(ps_18S_raw_2_percent)[,"Division"], exclude = NULL) #check that Eukaryota_X has been removed

#9.5. REMOVE UNWANTED SAMPLES ####
#ASV count per sample 
readcount(ps_18S_raw_2_percent)
#This shows that BYL48.2, MKIS3, CBP3C, and CBP9B all contain very low count data and are therefore removed for downstream analysis

#NOTE THIS IS AN UPDATED SCRIPT YOU MADE IN 30/03/23 to see what difference prevalence filtering makes on the data
ps_18S_new_trim <- ps_18S_raw_2_percent #now only 5664 ASVs

ps_18S_new_trim <- prune_samples(sample_names(ps_18S_new_trim) != "BYL48.2", ps_18S_new_trim)
ps_18S_new_trim <- prune_samples(sample_names(ps_18S_new_trim) != "MKIS3", ps_18S_new_trim)
ps_18S_new_trim <- prune_samples(sample_names(ps_18S_new_trim) != "CBP3C", ps_18S_new_trim)
ps_18S_new_trim <- prune_samples(sample_names(ps_18S_new_trim) != "CBP9B", ps_18S_new_trim)
ps_18S_new_trim #should now have four less samples (93)

ps_18S_new_trim #5664 ASVs and 93 samples



#THIS IS WHERE I STOPPED ON 30/03/23 I DIDN'T RESAVE THE PS_18s_NEW_TRIM FILES JUST IN CASE...




#remove your DNA files again (you've already saved your raw data backup above)
dna_18S_new_trim <- Biostrings::DNAStringSet(taxa_names(ps_18S_new_trim))
names(dna_18S_new_trim) <- taxa_names(ps_18S_new_trim)
ps_18S_new_trim <- merge_phyloseq(ps_18S_new_trim, dna_18S_new_trim)
taxa_names(ps_18S_new_trim) <- paste0("ASV", seq(ntaxa(ps_18S_new_trim)))
ps_18S_new_trim
#Access this data:
refseq(ps_18S_new_trim)

#9.5. SAVE YOUR FILTERED PHYLOSEQ OBJECT
saveRDS(dna_18S_new_trim, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/filt_dna_18S_new_trim")
saveRDS(ps_18S_new_trim, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/filt_ps_18S_new_trim")

ps_18S_new_trim = readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/filt_ps_18S_new_trim")
ps_18S_new_trim
#9.6. EXPORT A FINAL SEQ AND TAX TABLE FROM PHYLOSEQ FOR USE IN EXCEL CONTAINING METAZOA AND STREPTOPHYTES#### 
#export an ASV count table
ps_asv_count_18S <- abundances(ps_18S_new_trim)
ps_asv_count_18S <- as.data.frame(ps_asv_count_18S)
ps_asv_count_18S <- ps_asv_count_18S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_count_18S, file = "ps_ASV_count_18S_new_trim.txt", sep = "\t", quote=F, col.names=NA)

#export a proportional abundance table of ASVs
ps_asv_ra_18S <- abundances(ps_18S_new_trim, "compositional")
ps_asv_ra_18S <- as.data.frame(ps_asv_ra_18S)
ps_asv_ra_18S <- ps_asv_ra_18S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_ra_18S, file = "ps_ASV_proportional_count_18S_new_trim.txt", sep = "\t", quote=F, col.names=NA)

#export a tax table
ps_tax <- tax_table(ps_18S_new_trim)
ps_tax <- as.data.frame(ps_tax)
ps_tax <- ps_tax %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_tax, file = "ps_18S_taxa_new_trim.txt", sep = "\t", quote=F, col.names=NA)

#export a merged ASV count and taxonomic assignment table
ASV_TAX_table_18S <- merge(ps_asv_count_18S, ps_tax, by.x = c("ASVID"), by.y = c("ASVID"))
write.table(ASV_TAX_table_18S, "merged_asv_tax_table_18S_new_trim.txt", sep = "\t", quote=, col.names=NA)

#export a merged relative abundance ASV count and taxonomic assignment table
ASV_TAX_table_18S_ra <- merge(ps_asv_ra_18S, ps_tax, by.x = c("ASVID"), by.y = c("ASVID"))
write.table(ASV_TAX_table_18S_ra, "merged_asv_tax_table_18S_new_trim_RA.txt", sep = "\t", quote=, col.names=NA)

#export your ASV sequences as a FASTA file
ASV_seqs <- refseq(ps_18S_new_trim)
head(ASV_seqs)
write.table(ASV_seqs, file = "18S_ASVs_new_trim.fasta", sep = "\t", quote=F, col.names=NA)

#export read count
read_count <- readcount(ps_18S_new_trim)

#Add this to your final phyloseq metadata object
sample_data(ps_18S_new_trim)$read_count <- read_count
sample_data(ps_18S_new_trim)

#psmelt function allows you to convert phyloseq objects into R dataframes, useful for future work with ggplot
meta <- sample_data(ps_18S_new_trim)
ps_meta <- meta(ps_18S_new_trim)
row.names(ps_meta)
#Write it to a .txt file
write.table(ps_meta, file = "ps_metadata_new_trim.txt", sep = "\t", quote=F, col.names=NA)

###I would recommend looking at your exported files on excel to understand their distribution between different samples###

#9.8. CREATE A PHYLOGENETIC TREE ####
#I struggled to run this locally, if you have a large ASV dataset (>1000 ASVs), I recommend running it on the command line using the script "DECIPHER_script.R"
ps_decipher_tree <- readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/ps_merged_18s_new_trim_2_4_22")

ps_decipher_tree 
plot_tree(ps_decipher_tree)

#9.9. PHYLOGENETIC AGGLOMERATION####
#This step uses the tip_glom feature of phyloseq to cluster ASVs based on their phylogenetic relatedness
#This also had to be carried out on a server using the script "glom_script.R"

#I tried this at 0.2, 0.4, 0.02, and 0.01 and selected 0.01 glom
h1 = 0.01

ps_glom <- readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/ps_merged_18S_glom_new_trim_0_01_4_4")

ps_glom
#output:
#otu_table()   OTU Table:         [ 9393 taxa and 93 samples ]
#sample_data() Sample Data:       [ 93 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 9393 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 9393 tips and 9391 internal nodes ]
#refseq()      DNAStringSet:      [ 9393 reference sequences ]

#Look at how agglomeration changes the phylogenetic tree
test_counts = plot_tree(ps_decipher_tree, method = "treeonly", ladderize = "left", title = "18S no agglomeration") + theme(plot.title = element_text(size = 15))
tipglom_counts_0.01 = plot_tree(ps_glom, method = "treeonly", ladderize = "left", title = "18S Glom 0.01") + theme(plot.title = element_text(size = 15))

grid.arrange(nrow = 1, test_counts, tipglom_counts_0.01)

#9.10. EXPORT YOUR FINAL PHYLOGENETIC AGGLOMERATED DATA TABLE
#Separate ref seqs
refseq(ps_glom)

#export an ASV count table
ps_asv_count_18S <- abundances(ps_glom)
ps_asv_count_18S <- as.data.frame(ps_asv_count_18S)
ps_asv_count_18S <- ps_asv_count_18S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_count_18S, file = "FINAL_ps_ASV_count_18S_glom_0_01.txt", sep = "\t", quote=F, col.names=NA)

#export a proportional abundance table of ASVs
ps_asv_ra_18S <- abundances(ps_glom,"compositional")
ps_asv_ra_18S <- as.data.frame(ps_asv_ra_18S)
ps_asv_ra_18S <- ps_asv_ra_18S %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_ra_18S, file = "FINAL_ps_ASV_proportional_count_18S_glom_0_01.txt", sep = "\t", quote=F, col.names=NA)

#export a tax table
ps_tax <- tax_table(ps_glom)
ps_tax <- as.data.frame(ps_tax)
ps_tax <- ps_tax %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_tax, file = "FINAL_ps_18S_taxa_glom_0_01.txt", sep = "\t", quote=F, col.names=NA)

#export a merged ASV count and taxonomic assignment table
ASV_TAX_table_18S <- merge(ps_asv_count_18S, ps_tax, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_18S
write.table(ASV_TAX_table_18S, "FINAL_merged_asv_tax_table_18S_glom_0_01.txt", sep = "\t", quote=, col.names=NA)

##export a merged relative abundance (proportional) ASV count and taxonomic assignment table
ASV_TAX_table_18S_ra <- merge(ps_asv_ra_18S, ps_tax, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_18S_ra
write.table(ASV_TAX_table_18S_ra, "FINAL_merged_asv_tax_table_18S_glom_0_01_RA.txt", sep = "\t", quote=, col.names=NA)

#export your ASV sequences as a FASTA file
ASV_seqs <- refseq(ps_glom)
ASV_seqs #check fasta correct number of sequences 
write.table(ASV_seqs, file = "FINAL_18S_ASVs_glom_0_01.fasta", sep = "\t", quote=F, col.names=NA)

#export read count
read_count <- readcount(ps_glom)

#Add the read count to your final phyloseq metadata object
sample_data(ps_glom)$read_count <- read_count
sample_data(ps_glom)

#psmelt function allows you to convert phyloseq objects into R dataframes, useful for future work with ggplot
meta <- sample_data(ps_glom)
write.table(meta, "FINAL_ps_metadata_18S_glom_0_01.txt", sep = "\t", quote=, col.names=NA)

#Save your final version of your phyloseq object for easy reloading later on
ps_glom
saveRDS(ps_glom, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FINAL_ps_18S")

#Create a histogram of your final ASV and Sample counts
readsumsdf = data.frame(nreads = sort(taxa_sums(ps_glom), TRUE), sorted = 1:ntaxa(ps_glom), type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps_glom), 
                                                        TRUE), sorted = 1:nsamples(ps_glom), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#Create a plot of the final prevalence of ASVs faceted taxonomically by Division 
prevdf_18S = apply(X = otu_table(ps_glom),
                       MARGIN = ifelse(taxa_are_rows(ps_glom), yes = 1, no = 2),
                       FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_18S = data.frame(Prevalence = prevdf_18S,
                            TotalAbundance = taxa_sums(ps_glom),
                            tax_table(ps_glom))
head(prevdf_18S)

write.table(prevdf_18S, "FINAL_prevalence_table_18S.txt")
saveRDS(prevdf_18S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FINAL_prevdf_18S_glom")

#Are Divisions comprised of mostly low-prevalence features? compute the (1) average (mean) and (2) total (sum) prevalences of the features in each phylum
plyr::ddply(prevdf_18S, "Division", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
prevdf1 <- subset(prevdf_18S, Division %in% get_taxa_unique(ps_glom, "Division"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps_glom),color=Division)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Division) + theme(legend.position="none")

#CREATE A METAZOA ONLY PHYLOSEQ OBJECT####
extractMet = c("Metazoa")
metazoa = subset_taxa(ps_glom, Division %in% extractMet)
metazoa #should be 1190 taxa in 93 samples
saveRDS(metazoa, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FINAL_ps_metazoa")

##REMOVE METAZOA AND EMBRYOPHYTES FOR DOWNSTREAM STATS####
filterMet = c("Metazoa") #change to your specific Supergroup(s)
ps_glom_no_metazoa = subset_taxa(ps_glom, !Division %in% filterMet)
ps_glom_no_metazoa #8203
filterLandPlants = c("Embryophyceae")
ps_glom_no_metazoa_plants = subset_taxa(ps_glom_no_metazoa, !Class %in% filterLandPlants)
ps_glom_no_metazoa_plants

#Confirm you want to not include metazoa and embryophytes
ps <- ps_glom_no_metazoa_plants
ps # for ease when saving files below:

#EXPORT COUNT, RA COUNT, AND TAXA TABLE FOR DOWNSTREAM STATS WITHOUT METAZOA AND EMBRYOPHYTES

#Separate ref seqs
refseq(ps)

#export an ASV count table
ps_asv_count_18S_no_metazoa_plants <- abundances(ps)
ps_asv_count_18S_no_metazoa_plants <- as.data.frame(ps_asv_count_18S_no_metazoa_plants)
ps_asv_count_18S_no_metazoa_plants <- ps_asv_count_18S_no_metazoa_plants %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_count_18S_no_metazoa_plants, file = "FINAL_ps_ASV_count_18S_glom_0_01_no_met_plant.txt", sep = "\t", quote=F, col.names=NA)

#export a proportional abundance table of ASVs
ps_asv_ra_18S_no_metazoa_plants <- abundances(ps,"compositional")
ps_asv_ra_18S_no_metazoa_plants <- as.data.frame(ps_asv_ra_18S_no_metazoa_plants)
ps_asv_ra_18S_no_metazoa_plants <- ps_asv_ra_18S_no_metazoa_plants %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_ra_18S_no_metazoa_plants, file = "FINAL_ps_ASV_proportional_count_18S_glom_0_01_no_met_plant.txt", sep = "\t", quote=F, col.names=NA)

#export a tax table
ps_tax_no_metazoa_plants <- tax_table(ps)
ps_tax_no_metazoa_plants <- as.data.frame(ps_tax_no_metazoa_plants)
ps_tax_no_metazoa_plants <- ps_tax_no_metazoa_plants %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_tax_no_metazoa_plants, file = "FINAL_ps_18S_taxa_glom_0_01_no_met_plants.txt", sep = "\t", quote=F, col.names=NA)

#export a merged ASV count and taxonomic assignment table
ASV_TAX_table_18S_no_metazoa_plants <- merge(ps_asv_count_18S_no_metazoa_plants, ps_tax_no_metazoa_plants, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_18S_no_metazoa_plants
write.table(ASV_TAX_table_18S_no_metazoa_plants, "FINAL_merged_asv_tax_table_18S_glom_0_01_no_metazoa_plants.txt", sep = "\t", quote=, col.names=NA)

##export a merged relative abundance (proportional) ASV count and taxonomic assignment table
ASV_TAX_table_18S_no_metazoa_plants_ra <- merge(ps_asv_ra_18S_no_metazoa_plants, ps_tax_no_metazoa_plants, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_18S_no_metazoa_plants_ra
write.table(ASV_TAX_table_18S_no_metazoa_plants_ra, "FINAL_merged_asv_tax_table_18S_glom_0_01_no_metazoa_plants_RA.txt", sep = "\t", quote=, col.names=NA)

#export your ASV sequences as a FASTA file
ASV_seqs <- refseq(ps)
ASV_seqs #check fasta correct number of sequences 
write.table(ASV_seqs, file = "FINAL_18S_ASVs_glom_0_01_no_met_plants.fasta", sep = "\t", quote=F, col.names=NA)

#export read count
read_count <- readcount(ps)

#Add the read count to your final phyloseq metadata object
sample_data(ps)$read_count <- read_count
sample_data(ps)

#psmelt function allows you to convert phyloseq objects into R dataframes, useful for future work with ggplot
meta <- sample_data(ps)
write.table(meta, "FINAL_ps_metadata_18S_glom_0_01_no_met_plants.txt", sep = "\t", quote=, col.names=NA)

#Save your final version of your phyloseq object for easy reloading later on with no metazoa or plants
ps
saveRDS(ps, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FINAL_ps_18S_no_met_plants")

#Create a histogram of your final ASV and Sample counts
readsumsdf = data.frame(nreads = sort(taxa_sums(ps), TRUE), sorted = 1:ntaxa(ps), type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps), 
                                                        TRUE), sorted = 1:nsamples(ps), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#Create a plot of the final prevalence of ASVs faceted taxonomically by Division 
prevdf_18S_no_met_plants = apply(X = otu_table(ps),
                   MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_18S_no_met_plants = data.frame(Prevalence = prevdf_18S_no_met_plants,
                        TotalAbundance = taxa_sums(ps),
                        tax_table(ps))
head(prevdf_18S_no_met_plants)

plyr::ddply(prevdf_18S_no_met_plants, "Division", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
prevdf1_no_met_plant <- subset(prevdf_18S_no_met_plants, Division %in% get_taxa_unique(ps, "Division"))
ggplot(prevdf1_no_met_plant, aes(TotalAbundance, Prevalence / nsamples(ps),color=Division)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Division) + theme(legend.position="none")

#Save these final versions of your prevalence data without metazoa and plants
write.table(prevdf_18S_no_met_plants, "FINAL_prevalence_table_18S_no_metazoa_plants.txt")
saveRDS(prevdf_18S_no_met_plants, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FINAL_prevdf_18S_glom_no_met_plants")

#END OF SCRIPT####