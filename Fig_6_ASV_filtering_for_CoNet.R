#CoNet input table filtering 
#Author: Patrick M. Hooper
#Created: 14.11.22

#Background
#CoNet requires low abundance and prevalence ASVs to be removed from datatables, this script extracts a filtered abundance table from my 16S and 18S phyloseq objects

#Load packages
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(ggplot2); packageVersion("ggplot2")
library(knitr); packageVersion("knitr")
library(gridExtra); packageVersion("gridExtra")
library(microbiome); packageVersion("microbiome")
library(microbiomeutilities); packageVersion("microbiomeutilities")
library(Biostrings); packageVersion("Biostrings")
library(tibble); packageVersion("tibble")

#Load the averaged sample phyloseq objects, replicates from each pond are averaged across the three samples.
#The data is not currently relative abundance transformed

#Load 16S phyloseq####
ps_16S <- readRDS(file = "~/ps_16S_ave")
ps_16S #15250 taxa
sample_names(ps_16S)

#Load 18S phyloseq####
ps_18S <- readRDS(file = "~/ps_func_18S_ave")
ps_18S #9354 taxa
sample_names(ps_18S)

#Make sure the two phyloseq objects have samples ordered in the same way
#This uses the microViz function to order your ps_18S file by your ps_18S file
library(microViz); packageVersion("microViz")
sample_order <- sample_names(ps_16S)
ps_reorder(ps_18S, sample_order)

#check order
sample_names(ps_18S) == sample_names(ps_16S)

#Relative abundance transformed
ps_16S_ra <- transform_sample_counts(ps_16S, function(x){x / sum(x)})
ps_16S_ra

ps_18S_ra <- transform_sample_counts(ps_18S, function(x){x / sum(x)})
ps_18S_ra

#DATASET 1: 18S ONLY
#The total 18S dataset including the metazoa and plant groups included

#1. Abundance Filtering
#This is from the phyloseq pre-processing pages, this function keeps only those ASVs that have a mean abundance of greater than 1x10^-5 = 0.001% mean abundance
ps_18S_ra_a <- filter_taxa(ps_18S_ra, function(x) mean(x) > 1e-5, TRUE)
ps_18S_ra_a #4054 taxa

#2. Prevalence Filtering
#This is a handy feature to filter based on sample prevalence where the threshold is the number of samples that need to contain each ASV in a phyloseq object, e.g 0.25 = 25% of all samples, i.e. 8/32 of our samples
#remotes::install_github("vmikk/metagMisc")
library(metagMisc); packageVersion("metagMisc")

ps_18S_ra_a_p <- phyloseq_filter_prevalence(ps_18S_ra_a, prev.trh = 0.25, abund.trh = NULL) 
ps_18S_ra_a_p #539 taxa

#Let's also make a genus-agglomerated version while we're here
ps_18S_ra_a_p_genus <- tax_glom(ps_18S_ra_a_p, taxrank = "Genus")
ps_18S_ra_a_p_genus #183 taxa

#We also need to remove the function I added to this file, this allows us to add the best hit to the ASV names which stops us getting shared ASV names between our 16S and 18S files
tax_table(ps_18S_ra_a_p) <- tax_table(ps_18S_ra_a_p)[,1:7]
colnames(tax_table(ps_18S_ra_a_p)) #we now have 7 levels of taxonomic annotation
#Remember to note that best hit will change the name of our taxonomic levels from how PR2 defines them, but just remember the order is the same

#Now we add our best hit
ps_18S_ra_a_p <- microbiomeutilities::format_to_besthit(ps_18S_ra_a_p)
head(otu_table(ps_18S_ra_a_p))
colnames(tax_table(ps_18S_ra_a_p)) #Note how this hs changed it to "Domain"   "Phylum"   "Class"    "Order"    "Family"   "Genus"    "Species"  "best_hit

#Repeat for our genus agglomerated ASV table
tax_table(ps_18S_ra_a_p_genus) <- tax_table(ps_18S_ra_a_p_genus)[,1:7]
colnames(tax_table(ps_18S_ra_a_p_genus))
ps_18S_ra_a_p_genus <- microbiomeutilities::format_to_besthit(ps_18S_ra_a_p_genus)
head(otu_table(ps_18S_ra_a_p_genus))
colnames(tax_table(ps_18S_ra_a_p_genus))

#Finally, let's simplify the names of our two filtered phyloseq objects
#This is the abundance and prevalence filtered 18S ASVs with metazoa, fungi, rhodophyta, and embryophyta included
net_18S_ASVs <- ps_18S_ra_a_p
net_18S_ASVs
#otu_table()   OTU Table:         [ 539 taxa and 32 samples ]
#sample_data() Sample Data:       [ 32 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 539 taxa by 8 taxonomic ranks ]

net_18S_ASVs_genus <- ps_18S_ra_a_p_genus
net_18S_ASVs_genus
#otu_table()   OTU Table:         [ 183 taxa and 32 samples ]
#sample_data() Sample Data:       [ 32 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 183 taxa by 10 taxonomic ranks ]
#This doesn't contain BEST HIT

#DATASET 2: 18S ONLY WITHOUT METAZOA, EMBRYOPHYTES, AND RHODOPHYTES
net_protist_fungi_ASVs <- subset_taxa(net_18S_ASVs, !Class %in% c("Metazoa", "Rhodophyta") & !Order  %in% c("Embryophyceae"))
net_protist_fungi_ASVs <- transform_sample_counts(net_protist_fungi_ASVs, function(x){x / sum(x)})
net_protist_fungi_ASVs #we are redoing the relative abundance after removing the ASVs
#otu_table()   OTU Table:         [ 433 taxa and 32 samples ]
#sample_data() Sample Data:       [ 32 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 389 taxa by 8 taxonomic ranks ]

#Let's also make a genus-agglomerated version without metazoa and plants while we're here
net_protist_fungi_ASVs_genus <- tax_glom(net_protist_fungi_ASVs, taxrank = "Genus")
net_protist_fungi_ASVs_genus #167 taxa

#DATASET 3: 16S ONLY
ps_16S_ra_a <- filter_taxa(ps_16S_ra, function(x) mean(x) > 1e-5, TRUE)
ps_16S_ra_a_p <- phyloseq_filter_prevalence(ps_16S_ra_a, prev.trh = 0.25, abund.trh = NULL) 
ps_16S_ra_a_p #1859 taxa

#Let's also make a genus-agglomerated version while we're here
ps_16S_ra_a_p_genus <- tax_glom(ps_16S_ra_a_p, taxrank = "Genus")
ps_16S_ra_a_p_genus #289 taxa

#Best Hit application
tax_table(ps_16S_ra_a_p) <- tax_table(ps_16S_ra_a_p)[,1:7]
colnames(tax_table(ps_16S_ra_a_p)) #we now have 7 levels of taxonomic annotation, "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

#Now we add our best hit
ps_16S_ra_a_p <- microbiomeutilities::format_to_besthit(ps_16S_ra_a_p)
head(otu_table(ps_16S_ra_a_p))

#Repeat Best Hit for the genus phyloseq object
tax_table(ps_16S_ra_a_p_genus) <- tax_table(ps_16S_ra_a_p_genus)[,1:7]
colnames(tax_table(ps_16S_ra_a_p_genus))
ps_16S_ra_a_p_genus <- microbiomeutilities::format_to_besthit(ps_16S_ra_a_p_genus)
head(otu_table(ps_16S_ra_a_p_genus))
ps_16S_ra_a_p_genus #289 taxa

net_16S_ASVs <- ps_16S_ra_a_p
net_16S_ASVs_genus <-ps_16S_ra_a_p_genus

#Now we need to export this data! 
#CoNet requires a joined ASV and taxonomy table, we have to make a few edits to it

#net_18S_ASVs
net_18S_count <- as.data.frame(abundances(net_18S_ASVs))
net_18S_count <- net_18S_count %>% rownames_to_column(var = "ASVID")
net_18S_TAX <- as.data.frame(tax_table(net_18S_ASVs))
net_18S_TAX <- net_18S_TAX %>% rownames_to_column(var = "ASVID")
net_18S_input_table <- merge(net_18S_count, net_18S_TAX, by.x = c("ASVID"), by.y = c("ASVID"))
write.table(net_18S_input_table, "net_18S_input_table.tsv", sep = "\t", quote=, col.names=NA)

#net_18S_ASVs_genus
net_18S_genus_count <- as.data.frame(abundances(net_18S_ASVs_genus))
net_18S_genus_count <- net_18S_genus_count %>% rownames_to_column(var = "ASVID")
net_18S_genus_TAX <- as.data.frame(tax_table(net_18S_ASVs_genus))
net_18S_genus_TAX <- net_18S_genus_TAX %>% rownames_to_column(var = "ASVID")
net_18S_genus_input_table <- merge(net_18S_genus_count, net_18S_genus_TAX, by.x = c("ASVID"), by.y = c("ASVID"))
write.table(net_18S_genus_input_table, "net_18S_genus_input_table.tsv", sep = "\t", quote=, col.names=NA)

#net_protist_fungi_ASVs
net_protist_fungi_count <- as.data.frame(abundances(net_protist_fungi_ASVs))
net_protist_fungi_count <- net_protist_fungi_count %>% rownames_to_column(var = "ASVID")
net_protist_fungi_TAX <- as.data.frame(tax_table(net_protist_fungi_ASVs))
net_protist_fungi_TAX <- net_protist_fungi_TAX %>% rownames_to_column(var = "ASVID")
net_protist_fungi_input_table <- merge(net_protist_fungi_count, net_protist_fungi_TAX, by.x = c("ASVID"), by.y = c("ASVID"))
write.table(net_protist_fungi_input_table, "net_protist_fungi_input_table.tsv", sep = "\t", quote=, col.names=NA)

#net_protist_fungi_ASVs_genus
net_protist_fungi_genus_count <- as.data.frame(abundances(net_protist_fungi_ASVs_genus))
net_protist_fungi_genus_count <- net_protist_fungi_genus_count %>% rownames_to_column(var = "ASVID")
net_protist_fungi_genus_TAX <- as.data.frame(tax_table(net_protist_fungi_ASVs_genus))
net_protist_fungi_genus_TAX <- net_protist_fungi_genus_TAX %>% rownames_to_column(var = "ASVID")
net_protist_fungi_genus_input_table <- merge(net_protist_fungi_genus_count, net_protist_fungi_genus_TAX , by.x = c("ASVID"), by.y = c("ASVID"))
write.table(net_protist_fungi_genus_input_table, "net_protist_fungi_genus_input_table.tsv", sep = "\t", quote=, col.names=NA)

#net_16S_ASVs
net_16S_count <- as.data.frame(abundances(net_16S_ASVs))
net_16S_count <- net_16S_count %>% rownames_to_column(var = "ASVID")
net_16S_TAX <- as.data.frame(tax_table(net_16S_ASVs))
net_16S_TAX <- net_16S_TAX %>% rownames_to_column(var = "ASVID")
net_16S_input_table <- merge(net_16S_count, net_16S_TAX, by.x = c("ASVID"), by.y = c("ASVID"))
write.table(net_16S_input_table, "net_16S_input_table.tsv", sep = "\t", quote=, col.names=NA)

#net_16S_ASVs_genus
net_16S_genus_count <- as.data.frame(abundances(net_16S_ASVs_genus))
net_16S_genus_count <- net_16S_genus_count %>% rownames_to_column(var = "ASVID")
net_16S_genus_TAX <- as.data.frame(tax_table(net_16S_ASVs_genus))
net_16S_genus_TAX <- net_16S_genus_TAX %>% rownames_to_column(var = "ASVID")
net_16S_genus_input_table <- merge(net_16S_genus_count, net_16S_genus_TAX, by.x = c("ASVID"), by.y = c("ASVID"))
write.table(net_16S_genus_input_table, "net_16S_genus_input_table.tsv", sep = "\t", quote=, col.names=NA)

##Format for CoNet in Excel##
#There was a lot of changes made to the Excel to make it formatted for CoNet, these are provided below. However, the edited versions are also provided in the 'data_files' folders.
                           
#The first column (containing an ascending numerical list) had to be removed from the joined ASV and tax tables for both the 16S and 18S in Excel
#I also changed the taxonomic annotation so it was one column seperated by ";"s, using the concat rows feature in Excel
#4. Change the ASVID column header on column 1 to #OTUID
#5. Change the colons introduced by Best HIT to the ASVIDs column 1 to _ using find and replace
#6. Change the 'Domain' header to ''taxonomy' and remove all other taxonomic headings
#7. Delete the final column with the best_hit ID
#8. #I then checked any that already had a prefix that needed to be deleted, this takes a while! (find and replace '__D' with '')
#9. FINALLY, This is a bit more complicated, we need to add the SILVA style annotation to our taxonomy:
#k__KINGDOM;p__PHYLUM;c__CLASS;o__ORDER;f__FAMILY;g__GENUS;s__SPECIES,
#I did this with concatenate in excel by making a column for each prefix and then merging all the rows
#e.g. k__	Bacteria	p__	Firmicutes	c__	Clostridia	o__	Clostridiales	f__	Clostridiaceae	g__	Clostridium_sensu_stricto_13	s__	g__Clostridium_sensu_stricto_13
#Then concatenate it to make one column with all this information in, then you're done!
#Also remove any irregular characters like brackets, sometimes things have these, e.g. SAR3324_clade(Marine_group_B)
#SAVED THE FINAL VERSIONS OF THIS IN THE NETWORK PAGE AS AN EXCEL AND A TAB-DELIMITED FILE

####End of Script####
