#Assembly of V9 18S sequence data using dada2
#Author: Patrick M. Hooper
#Created: 01.04.22

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

##YOU'RE READY TO BEGIN WITH DADA2
#My sequencing files are from two separate sequencing runs, DADA2 requires you to assemble sequencing runs separately and merge them after assembly
#I have two sequencing runs from 2016 and 2018, I started with the 2016 sequences

#1.1. IDENTIFY PATH TO YOUR 18S FASTA SEQUENCING FILES####
path <- "~/2016_18S_Files"
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

#2.1. FILTER AND TRIM####
# Place filtered files in 'filtered/' subdirectory
filtFs_18S_2016 <- file.path(path, "filtered", paste0(sample.names_18S_2016, "_F_filt.fastq.gz"))
filtRs_18S_2016 <- file.path(path, "filtered", paste0(sample.names_18S_2016, "_R_filt.fastq.gz"))
names(filtFs_18S_2016) <- sample.names_18S_2016
names(filtRs_18S_2016) <- sample.names_18S_2016

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

#3.1. LEARN ERROR RATES####
errF_18S_2016 <- learnErrors(filtFs_18S_2016, multithread = TRUE)
errR_18S_2016 <- learnErrors(filtRs_18S_2016, multithread = TRUE)
errF_18S_2016
errR_18S_2016

#3.2. Plot your error graphs
#The error rates for each possible nucleotide transition are shown. Points = error rates for each quality score, black line = estimated error rates from the machine-learning algorithm. Red lines shows expected error rate under nominal Q-score.
#INTERPRETATION: If points fit well to the black line and error rate drops with increased quality score you can assume a good error rate.

Plot_Error_FWD <- plotErrors(errF_18S_2016, nominalQ=TRUE)#FWD read errors
Plot_Error_REV <- plotErrors(errR_18S_2016, nominalQ=TRUE) #REV read errors
Plot_Error_FWD
Plot_Error_REV

#4.1. SAMPLE INFERENCE####
#Apply the core sample inference algorithm to the filtered and trimmed sequence data
#Info: https://www.nature.com/articles/nmeth.3869#methods
dadaFs_18S_2016 <- dada(filtFs_18S_2016, err = errF_18S_2016, multithread = FALSE)
dadaRs_18S_2016 <- dada(filtRs_18S_2016, err = errR_18S_2016, multithread = FALSE)

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

#5.1. CONSTRUCT SEQUENCE TABLE####
seqtab_18S_2016 <- makeSequenceTable(mergers_18S_2016)
dim(seqtab_18S_2016)

#5.2. Inspect distribution of sequence lengths
seq_length_18S_2016 = table(nchar(getSequences(seqtab_18S_2016)))
seq_length_18S_2016
write.csv(seq_length_18S_2016, file = "~/18S_2016/seq_length_distribution")

#5.3 Trim those sequences that are significantly larger than the expected primer length ~130 for V9 
seqtab_18S_2016_trim<- seqtab_18S_2016[,nchar(colnames(seqtab_18S_2016)) %in% seq(95,150)]
table(nchar(getSequences(seqtab_18S_2016_trim)))

#5.3. REMOVE CHIMERAS
seqtab_nochim_18S_2016 <- removeBimeraDenovo(seqtab_18S_2016_trim, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab_nochim_18S_2016)
table(nchar(getSequences(seqtab_nochim_18S_2016)))

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

#Outside of filtering, there should no step in which a majority of reads are lost
#If a majority of reads failed to merge, you may need to revisit the truncLen parameter used in the filtering step and make sure that the truncated reads span your amplicon
#If a majority of reads were removed as chimeric, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification

#OPTIONAL 2: IDENTIFY PATH TO YOUR SECOND FOLDER OF SEQUENCING FILES####
#Skip if you're only analysing one sequencing run
path <- "~/2018_18S_Files"
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

#SECTION 3: Learn the error rates
errF_18S_2018 <- learnErrors(filtFs_18S_2018, multithread = TRUE)
errR_18S_2018 <- learnErrors(filtRs_18S_2018, multithread = TRUE)
errF_18S_2018
errR_18S_2018

#The error rates for each possible nucleotide transition are shown. Points = error rates for each quality score, black line = estimated error rates from the machine-learning algorithm. Red lines shows expected error rate under nominal Q-score. If points fit well to the black line and error rate drops with increased quality score you can assume a good error rate.

Plot_Error_FWD_2018 <- plotErrors(errF_18S_2018, nominalQ=TRUE) #FWD read errors
Plot_Error_REV_2018 <- plotErrors(errR_18S_2018, nominalQ=TRUE) #REV read errors
Plot_Error_FWD_2018
Plot_Error_REV_2018

#SECTION 4: Sample inference
#Apply the core sample inference algorithm to the filtered and trimmed sequence data
#Info: https://www.nature.com/articles/nmeth.3869#methods
dadaFs_18S_2018 <- dada(filtFs_18S_2018, err=errF_18S_2018, multithread=FALSE)
dadaRs_18S_2018 <- dada(filtRs_18S_2018, err=errR_18S_2018, multithread=FALSE)

dadaFs_18S_2018[[1]]
dadaRs_18S_2018[[1]]

#SECTION 5: MERGE YOUR SAMPLE INFERENCE WITH YOUR SEQUENCE READS
mergers_18S_2018 <- mergePairs(dadaFs_18S_2018, filtFs_18S_2018, dadaRs_18S_2018, filtRs_18S_2018, verbose=TRUE)
head(mergers_18S_2018)

#Inspect the merger data.frame from the first sample
head(mergers_18S_2018[[1]])

#SECTION 6: CONSTRUCT SEQUENCE TABLE
seqtab_18S_2018 <- makeSequenceTable(mergers_18S_2018)
dim(seqtab_18S_2018)

#5.2. Inspect distribution of sequence lengths
seq_length_18S_2018 = table(nchar(getSequences(seqtab_18S_2018)))
seq_length_18S_2018
write.csv(seq_length_18S_2018, file = "~/18S_2018/seq_length_distribution")

#5.3 Trim those sequences that are significantly larger than the expected primer length ~130 for V9 
seqtab_18S_2018_trim<- seqtab_18S_2018[,nchar(colnames(seqtab_18S_2018)) %in% seq(95,150)]
table(nchar(getSequences(seqtab_18S_2018_trim)))

#check sample names
rownames(seqtab_18S_2018)

#REMOVE CHIMERAS
seqtab_nochim_18S_2018 <- removeBimeraDenovo(seqtab_18S_2018_trim, method="consensus", multithread=FALSE, verbose=TRUE)

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

write.table(track_18S_2018, file = "track_18S_2018_new_trim.txt")

#IF YOU WERE WORKING ON ONE SEQUENCING RUN START AGAIN HERE#

#7.1. COMBINE YOUR SEQUENCING RUNS ####
#Check you've got ALL your samples and no replicate names
seqtab_nochim_18S_2016 = readRDS(file = "~/18S_2016/seqtab_nochim_18S_2016")
row.names(seqtab_nochim_18S_2016)
row.names(seqtab_nochim_18S_2018) 

#IF YOU HAVE MULTIPLE SEQUENCE RUNS THIS IS THE TIME YOU MERGE THEM TOGETHER BEFORE ASSIGNING TAXONOMY, IF YOU'RE JUST WORKING WITH ONE SEQUENCE RUN IGNORE THIS STEP
st_18S <- mergeSequenceTables(seqtab_nochim_18S_2016, seqtab_nochim_18S_2018)
head(st_18S)
row.names(st_18S)

# Inspect distribution of sequence lengths
table(nchar(getSequences(st_18S)))

#8.1. ASSIGN TAXONOMY WITH PR2 ####
merged_taxa_18S <- assignTaxonomy(st_18S, "pr2_version_4.14.0_SSU_dada2.fasta.gz", multithread=TRUE, taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"))

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
#Make a phyloseq object from your sequence table, tax file, and metadata.
ps_18S <- phyloseq(
  otu_table(st_18S, taxa_are_rows = FALSE),
  #taxa_are_rows is FALSE if importing directly from dada2 but TRUE if importing from QIIME2
  sample_data(metadata_18S),
  tax_table(merged_taxa_18S)
)
ps_18S #check your phyloseq object, how many ASVs and samples
sample_names(ps_18S) 

#9.2. FILTERING ASV TABLES IN PHYLOSEQ ####
#Investigate the total reads per sample and their distribution using a histogram
readsumsdf = data.frame(nreads = sort(taxa_sums(ps_18S), TRUE), sorted = 1:ntaxa(ps_18S_new_trim), 
                        type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps_18S), 
                                                        TRUE), sorted = 1:nsamples(ps_18S), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#9.4. TAXONOMIC FILTERING ####
rank_names(ps_18S) # show available rank names in the dataset

# Create table, number of features for each taxonomic division
table(tax_table(ps_18S)[, "Division"], exclude = NULL)

#remove any bacteria ASVs
ps_18S = subset_taxa(ps_18S, Kingdom == "Eukaryota")
ps_18S

#remove anything that is considered ambiguous annotation, i.e. NA annotation at Supergroup level
ps_18S <- subset_taxa(ps_18S, !is.na(Supergroup) & !Supergroup %in% c("Eukaryota_X", "NA"))
ps_18S

#create a table of the number of ASVs in each Division (L2)
table(tax_table(ps_18S)[,"Division"], exclude = NULL) 
#This shows us which divisions only feature a few ASVS, these might be worth filtering

#Define phyla to filter such as anything that appears in a very low percent of all samples or plastid and mito
filterDivision = c("Stramenopiles_X", "Picozoa", "Telonemia") #change to your specific Supergroup(s)

# Remove filtered Divisions
ps_18S = subset_taxa(ps_18S, !Division %in% filterDivision)

table(tax_table(ps_18S)[,"Division"], exclude = NULL) #check that Eukaryota_X has been removed

#9.5. REMOVE UNWANTED SAMPLES ####
#ASV count per sample 
readcount(ps_18S)
#This shows that BYL48.2, MKIS3, CBP3C, and CBP9B all contain very low count data and are therefore removed for downstream analysis

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

#CREATE A METAZOA ONLY PHYLOSEQ OBJECT####
extractMet = c("Metazoa")
metazoa = subset_taxa(ps_glom, Division %in% extractMet)
metazoa #should be 1190 taxa in 93 samples
saveRDS(metazoa, file = "~/FINAL_ps_metazoa")

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

#export an ASV count table
ps_asv_count_18S_no_metazoa_plants <- abundances(ps)
ps_asv_count_18S_no_metazoa_plants <- as.data.frame(ps_asv_count_18S_no_metazoa_plants)
ps_asv_count_18S_no_metazoa_plants <- ps_asv_count_18S_no_metazoa_plants %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_count_18S_no_metazoa_plants, file = "ps_18S_no_met_plant.txt", sep = "\t", quote=F, col.names=NA)

#export a proportional abundance table of ASVs
ps_asv_ra_18S_no_metazoa_plants <- abundances(ps,"compositional")
ps_asv_ra_18S_no_metazoa_plants <- as.data.frame(ps_asv_ra_18S_no_metazoa_plants)
ps_asv_ra_18S_no_metazoa_plants <- ps_asv_ra_18S_no_metazoa_plants %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_asv_ra_18S_no_metazoa_plants, file = "ps_proportional_count_18S_no_met_plant.txt", sep = "\t", quote=F, col.names=NA)

#export a tax table
ps_tax_no_metazoa_plants <- tax_table(ps)
ps_tax_no_metazoa_plants <- as.data.frame(ps_tax_no_metazoa_plants)
ps_tax_no_metazoa_plants <- ps_tax_no_metazoa_plants %>% 
  rownames_to_column(var = "ASVID")
write.table(ps_tax_no_metazoa_plants, file = "ps_18S_no_met_plants.txt", sep = "\t", quote=F, col.names=NA)

#export a merged ASV count and taxonomic assignment table
ASV_TAX_table_18S_no_metazoa_plants <- merge(ps_asv_count_18S_no_metazoa_plants, ps_tax_no_metazoa_plants, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_18S_no_metazoa_plants
write.table(ASV_TAX_table_18S_no_metazoa_plants, "merged_asv_tax_table_18S_no_metazoa_plants.txt", sep = "\t", quote=, col.names=NA)

##export a merged relative abundance (proportional) ASV count and taxonomic assignment table
ASV_TAX_table_18S_no_metazoa_plants_ra <- merge(ps_asv_ra_18S_no_metazoa_plants, ps_tax_no_metazoa_plants, by.x = c("ASVID"), by.y = c("ASVID"))
ASV_TAX_table_18S_no_metazoa_plants_ra
write.table(ASV_TAX_table_18S_no_metazoa_plants_ra, "merged_asv_tax_table_18S_no_metazoa_plants_RA.txt", sep = "\t", quote=, col.names=NA)

#export your ASV sequences as a FASTA file
ASV_seqs <- refseq(ps)
ASV_seqs #check fasta correct number of sequences 
write.table(ASV_seqs, file = "18S_ASVs_no_met_plants.fasta", sep = "\t", quote=F, col.names=NA)

#export read count
read_count <- readcount(ps)

#Add the read count to your final phyloseq metadata object
sample_data(ps)$read_count <- read_count
sample_data(ps)

#psmelt function allows you to convert phyloseq objects into R dataframes, useful for future work with ggplot
meta <- sample_data(ps)
write.table(meta, "ps_metadata_18S_no_met_plants.txt", sep = "\t", quote=, col.names=NA)

#Save your final version of your phyloseq object for easy reloading later on with no metazoa or plants
ps
saveRDS(ps, file = "~/FINAL_ps_18S_no_met_plants")

####END OF SCRIPT####
