#V4 16S and V9 18S ordination 
#Author: Patrick M. Hooper
#Created: 30/08/22
#Updated: 08/09/22

#1. LOAD PACKAGES####
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(Biostrings); packageVersion("Biostrings")
library(microbiome); packageVersion("microbiome")
library(DESeq2); packageVersion("DESeq2")
library(vegan); packageVersion("vegan")
library(microbiomeutilities); packageVersion("microbiomeutilities")

#OPTIONAL: set ggplot2 theme to minimal
theme_set(theme_minimal()) 

#Set standardized plot font sizes
font_size <- theme(axis.text = element_text(size = 18)) + theme(axis.title = element_text(size = 18)) + theme(plot.title = element_text(size = 20)) + theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size = 20)) + theme(plot.subtitle = element_text(size = 18))

##OPTIONAL: INSTALL AND LOAD PALETTES##
library("chroma")
library("viridis")

#2. SET WORKING DIRECTORY####
setwd("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS")

#3. Set palette####
colour_palette <- c("1_kuujjuarapik" = "#228833","2_umiujaq" = "#CCBB44", "3_cambridge_bay" = "#EE6677","4_bylot_island" = "#AA3377", "5_resolute" = "#66CCEE","6_ellesmere_island" = "#4477AA","7_ice_shelves" = "#BBBBBB") 
show_col(colour_palette)

#4. Loading 18S phyloseq objects
#We are using count tables averaged by water bodies. Biological replicates from water bodies were summed and averaged by the number of replicates. The metadata was also averaged. 

#5. Add averaged metadata
#I originally debated changing all the IDs on my file headings but that is a bad idea as its much easier to just add the 'new_IDs' as extra metadata in the sample_data spreadsheet rather than manually edit my sheets
#Generally going forward it's best to do as much as possible within R rather than in excel

#this table includes pond data to group the samples and the new sample IDs for the plots
#No KJB17bii in this file
meta_18S <- read.table("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/GLOM_V9_18S_METADATA_FORMATTED_UPDATE_AVERAGE_17_8_22.txt") #NOTE THIS IS THE METADATA WITH THE READ COUNTS WITHOUT METAZOA AND PLANTS
meta_18S
rownames(meta_18S)
nrow(meta_18S) #32

#6. LOAD YOUR ASV COUNT TABLE .TXT FILE
#I made this in excel from the functional mothership to make sure that ASVIDs were all matching
#No KJB17bii in this file
#This needs to not have a column header on the ASVID column
count_18S <- read.table("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/FINAL_ps_18S_asv_count_glom_0_01_FUNCTIONAL_AVERAGE_25_8.txt", header = TRUE) 

#I also had to change it to a matrix
count_18S <- as.matrix(count_18S)

#read column names
colnames(count_18S)
ncol(count_18S) #32

#read rownames
head(rownames(count_18S)) #check in numerical order

#total number of rows (ASVs)
nrow(count_18S) # 1 ASV was removed at this stage as it was not taxonomically annotated higher than eukaryotic

#number of reads in each sample
colSums(count_18S)

#7. LOAD YOUR FUNCTIONAL ASV TAXA TABLE WITH YOUR ADDED FUNCTIONAL DATA
#THIS WAS CREATED USING VLOOKUP AGAINST THE MAIN ANNOTATED FILE
#I had to manually remove the heading ASVID on the ASVs for phyloseq to judge it worthy for use
#No KJB17bii in this file

#I updated this following Alan Warren's recommendation in May 23
func_tax_18S <- read.table("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/18S_functional_taxonomy_table_08_12_22.txt", header = TRUE)
head(func_tax_18S)

#use this function to add ASVIDS as rownames on the tax table
rownames(func_tax_18S) <- func_tax_18S[,1]
rownames(func_tax_18S)

#Need to remove the column ASVID
func_tax_18S <- func_tax_18S %>% select(-ASVID)

#Check your function names:
unique(func_tax_18S$Function)

#check your nutrition columns:
unique(func_tax_18S$Nutrition)

#I also had to change it to a matrix
func_tax_18S <- as.matrix(func_tax_18S)

#number of rows in taxa table
nrow(func_tax_18S)

head(func_tax_18S) #this shows the rownames are now are ASV IDs

#Taxonomy levels:
colnames(func_tax_18S)
#Output: 
#[1] "Kingdom"    "Supergroup" "Division"   "Class"     
#[5] "Order"      "Family"     "Genus"      "Species"   
#[9] "Function"   "Func_Group"  "Nutrition" "Notes"

#Check all your ASVs match
class(count_18S)
class(func_tax_18S) # needs to be a matrix
all(rownames(count_18S) == rownames(func_tax_18S)) #TRUE = GOOD

#read in as phyloseq objects
OTU_18S_ave	=	otu_table(count_18S,	taxa_are_rows	=	TRUE)
head(OTU_18S_ave)
TAX_18S_ave	=	tax_table(func_tax_18S)
head(TAX_18S_ave)
META_18S_ave	=	sample_data(meta_18S)

#read in as phyloseq objects
ps_func_18S_ave = merge_phyloseq(OTU_18S_ave, TAX_18S_ave, META_18S_ave)
ps_func_18S_ave #You're good to go! You now have an averaged ASV count table of all your 18S reads, averaged metadata, and your functional info for your 18S reads

#remove 0 ASVs
#Check that there are no ASVs that are not present in any of the samples
ps_func_18S_ave <- prune_taxa(taxa_sums(ps_func_18S_ave) > 0, ps_func_18S_ave)
ps_func_18S_ave #this should leave us with 9354 ASVs by 92 samples

#phyloseq-Order experiment-level object
#otu_table()   OTU Table:         [ 9354 taxa and 32 samples ]
#sample_data() Sample Data:       [ 32 samples by 20 sample variables ]
#tax_table()   Taxonomy Table:    [ 9354 taxa by 9 taxonomic ranks ]

#quick check all names are correct
taxa_names(TAX_18S_ave)
taxa_names(OTU_18S_ave)
sample_names(META_18S_ave)
sample_names(OTU_18S_ave)

#save phyloseq object
saveRDS(ps_func_18S_ave, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/ps_func_18S_ave")

#4. Loading 16S phyloseq objects
#We are using count tables averaged by water bodies. Biological replicates from water bodies were summed and averaged by the number of replicates. The metadata was also averaged. 


#3. Add averaged metadata
#I originally debated changing all the IDs on my file headings but that is a bad idea as its much easier to just add the 'new_IDs' as extra metadata in the sample_data spreadsheet rather than manually edit my sheets
#Generally going forward it's best to do as much as possible within R rather than in excel

#this table includes pond data to group the samples and the new sample IDs for the plots
#No KJB17bii in this file
meta_16S <- read.table("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/GLOM_V4_16S_METADATA_FORMATTED_AVERAGE_UPDATE_17_8_22.txt") #NOTE THIS IS THE METADATA WITH THE READ COUNTS WITHOUT METAZOA AND PLANTS
head(meta_16S)
rownames(meta_16S)
nrow(meta_16S) #32

#6. LOAD YOUR ASV COUNT TABLE .TXT FILE
#I made this in excel from the functional mothership to make sure that ASVIDs were all matching
#No KJB17bii in this file
#This needs to not have a column header on the ASVID column
count_16S <- read.table("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/FINAL_ps_16S_count_glom_0_01_AVERAGE.txt", header = TRUE) 


#I also had to change it to a matrix
count_16S <- as.matrix(count_16S)

#read column names
colnames(count_16S)
ncol(count_16S) #32

#read rownames
head(rownames(count_16S)) #check in numerical order

#total number of rows (ASVs)
nrow(count_16S) #15250

#number of reads in each sample
colSums(count_16S)

#7. LOAD YOUR 16S ASV TAXA TABLE
#THIS WAS CREATED USING VLOOKUP AGAINST THE MAIN ANNOTATED FILE
#I had to manually remove the heading ASVID on the ASVs for phyloseq to judge it worthy for use
#No KJB17bii in this file
tax_16S <- read.table("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/FINAL_ps_16S_tax_no_KJ17bii.txt", header = TRUE) #had an issue here for a while as got the error "more columns than column names", this was due to spaces in the annotation given by SILVA database
head(tax_16S)


#use this function to add ASVIDS as rownames on the tax table
rownames(tax_16S) <- tax_16S[,1]
rownames(tax_16S)


#Need to remove the column ASVID
tax_16S <- tax_16S %>% select(-ASVID)


#I also had to change it to a matrix
tax_16S<- as.matrix(tax_16S)

#number of rows in taxa table
nrow(tax_16S)

head(tax_16S) #this shows the rownames are now are ASV IDs

#Taxonomy levels:
colnames(tax_16S)
#Output: 
#"Kingdom"    "Supergroup" "Division"   "Order"      "Order"     
#"Family"     "Genus"      "Species"   

#Check all your ASVs match
Order(count_16S)
Order(tax_16S) # needs to be a matrix
all(rownames(count_16S) == rownames(tax_16S)) #TRUE = GOOD

#read in as phyloseq objects
OTU_16S_ave	=	otu_table(count_16S,	taxa_are_rows	=	TRUE)
head(OTU)
TAX_16S_ave	=	tax_table(tax_16S)
head(TAX)
META_16S_ave	=	sample_data(meta_16S)

#read in as phyloseq objects
ps_16S_ave = merge_phyloseq(OTU_16S_ave, TAX_16S_ave, META_16S_ave)
ps_16S_ave #You're good to go! You now have an averaged ASV count table of all your 18S reads, averaged metadata

#remove 0 ASVs
#Check that there are no ASVs that are not present in any of the samples
ps_16S_ave <- prune_taxa(taxa_sums(ps_16S_ave) > 0, ps_16S_ave)
ps_16S_ave #this should leave us with 15250 by 32 samples

#phyloseq-Order experiment-level object
#otu_table()   OTU Table:         [ 15250 taxa and 32 samples ]
#sample_data() Sample Data:       [ 32 samples by 21 sample variables ]
#tax_table()   Taxonomy Table:    [ 15250 taxa by 7 taxonomic ranks ]

#quick check all names are correct
taxa_names(TAX_16S_ave)
taxa_names(OTU_16S_ave)
sample_names(META_16S_ave)
sample_names(OTU_16S_ave)

#save phyloseq object
saveRDS(ps_16S_ave, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/ps_16S_ave")

#transform sample counts to relative abundance (proportional) counts
ps_func_18S_ave_RA  = transform_sample_counts(ps_func_18S_ave, function(x) x / sum(x))
ps_16S_ave_RA = transform_sample_counts(ps_16S_ave, function(x) x / sum(x))
ps_func_18S_ave_RA
ps_16S_ave_RA

#Now you have loaded your data it is time to make some sub-divisions of your phyloseq for plotting

#Make protist, fungi and metazoa only tables
FUNGI_func_ave_RA = subset_taxa(ps_func_18S_ave_RA, Division == "Fungi")
METAZOA_func_ave_RA = subset_taxa(ps_func_18S_ave_RA, Division == "Metazoa")
PROTIST_func_ave_RA = subset_taxa(ps_func_18S_ave_RA, !Division %in% c("Fungi", "Metazoa","") & !Order  %in% c("Embryophyceae"))
FUNGI_func_ave_RA #1804 taxa ASVs, 32 samples
METAZOA_func_ave_RA #1184 taxa ASVs, 32 samples
PROTIST_func_ave_RA #6366 taxa ASVs

#Make cyanobacteria and proteobacteria only tables
cyanobacteria_ave_RA = subset_taxa(ps_16S_ave_RA, Phylum == "Cyanobacteria")
proteobacteria_ave_RA = subset_taxa(ps_16S_ave_RA, Phylum == "Proteobacteria")
cyanobacteria_ave_RA #734 taxa
proteobacteria_ave_RA #3658 taxa

#Bray-Curtis NMDS####
#1. 18S average
#2. protist average
#3. metazoa average
#4. fungi average
#5. 16S average
#6. cyanobacteria average
#7. proteobacteria average

#Why are we using NMDS?
#Taken from: https://jkzorz.github.io/2019/06/06/NMDS.html
#Non-metric Multi-dimensional Scaling (NMDS) is a way to condense information from multidimensional data (multiple variables/species/OTUs), into a 2D representation or ordination. In this ordination, the closer two points are, the more similar the corresponding samples are with respect to the variables that went into making the NMDS plot.

#NMDS plots are great tools for microbial ecologists (and others working with big data) because you can condense the overwhelming amount of information from the distribution of multiple species/OTUs across your samples into something you can actually look at. In an NMDS plot generated using an OTU table the points are samples. The closer the points/samples are together in the ordination space, the more similar their microbial communities.

#NMDS plots are non-metric, meaning that among other things, they use data that is not required to fit a normal distribution. This is handy for microbial ecologists because the majority of our data has a skewed distribution with a long tail. In other words, there are only a few abundant species, and many, many species with low abundance (the long tail).

#NMDS is an iterative algorithm, so it repeats the same series of steps over and over again until it finds the best solution. This is important to note because it means that each time you produce an NMDS plot from scratch it may look slightly different, even when starting with exactly the same data.

#What makes an NMDS plot non-metric is that it is rank-based. This means that instead of using the actual values to calculate distances, it uses ranks. So for example, instead of saying that sample A is 5 points away from sample B, and 10 from sample C, you would instead say that: sample A is the "1st" most close sample to B, and sample C is the "2nd" most close.

#Why Bray-Curtis 
#The last basic thing to know about NMDS is that it uses a distance matrix as an input. Read more about distance measures here. There are many different distance measures to choose from, however as a default, I tend to use Bray-Curtis when dealing with relative abundance data. Bray-Curtis takes into account species presence/absence, as well as abundance, whereas other measures (like Jaccard) only take into account presence/absence.
#It is invariant to changes in units.
#It is unaffected by additions/removals of species that are not present in two communities,
#It is unaffected by the addition of a new community,
#It can recognize differences in total abundances when relative abundances are the same

#Bray-Curtis NMDS ordinations
ps_func_18S_ave_RA_BCNMDS <- ordinate(ps_func_18S_ave_RA, distance = "bray", method = "NMDS")
ps_func_18S_ave_RA_BCNMDS
#Dimensions: 2 
#Stress:     0.1770384
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'veganifyOTU(physeq)

FUNGI_func_ave_RA_BCNMDS <- ordinate(FUNGI_func_ave_RA, distance = "bray", method = "NMDS")
FUNGI_func_ave_RA_BCNMDS
#Dimensions: 2 
#Stress:     0.2201285 
#Stress type 1, weak ties
#No convergent solutions - best solution after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'veganifyOTU(physeq)' 

METAZOA_func_ave_RA_BCNMDS <- ordinate(METAZOA_func_ave_RA, distance = "bray", method = "NMDS")
METAZOA_func_ave_RA_BCNMDS
#Dimensions: 2 
#Stress:     0.166802
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'veganifyOTU(physeq)' 

PROTIST_func_ave_RA_BCNMDS <- ordinate(PROTIST_func_ave_RA, distance = "bray", method = "NMDS")
PROTIST_func_ave_RA_BCNMDS
#Dimensions: 2 
#Stress:     0.1920405  
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'veganifyOTU(physeq)' 

ps_16S_ave_RA_BCNMDS <- ordinate(ps_16S_ave_RA, distance = "bray", method = "NMDS")
ps_16S_ave_RA_BCNMDS 
#Dimensions: 2 
#Stress:     0.1038522 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'veganifyOTU(physeq)' 

cyanobacteria_ave_RA_BCNMDS <- ordinate(cyanobacteria_ave_RA, distance = "bray", method = "NMDS")
cyanobacteria_ave_RA_BCNMDS
#Dimensions: 2 
#Stress:     0.1712805 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'veganifyOTU(physeq)' 

proteobacteria_ave_RA_BCNMDS <- ordinate(proteobacteria_ave_RA, distance = "bray", method = "NMDS")
proteobacteria_ave_RA_BCNMDS
#Dimensions: 2 
#Stress:    0.1046206  
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'veganifyOTU(physeq)' 

#Plotting your ordination - 18S ####
#Average 18S with functional annotation
p_ps_func_18S_ave_RA_BCNMDS <- plot_ordination(ps_func_18S_ave_RA, ps_func_18S_ave_RA_BCNMDS, type="samples", color="region", shape = "water_body") 
p_ps_func_18S_ave_RA_BCNMDS  <- p_ps_func_18S_ave_RA_BCNMDS + geom_point(size=8, alpha = 0.8) + ggtitle("Relative Abundance counts of 18S ASVs averaged by water body", subtitle = "Bray-Curtis NMDS Stress = 0.177, ANOSIM R-Stat 0.7046, p value < 0.001")
p_ps_func_18S_ave_RA_BCNMDS $layers
p_ps_func_18S_ave_RA_BCNMDS $layers <- p_ps_func_18S_ave_RA_BCNMDS $layers[-1]
p_ps_func_18S_ave_RA_BCNMDS <- p_ps_func_18S_ave_RA_BCNMDS + scale_colour_manual(values = colour_palette)
p_ps_func_18S_ave_RA_BCNMDS + font_size

#use ggsave to save to your wd as a pdf and jpeg
ggsave("p_ps_func_18S_ave_RA_BCNMDS.jpg", width = 297, height = 210, units = c("mm"))
ggsave("p_ps_func_18S_ave_RA_BCNMDS.pdf", width = 297, height = 210, units = c("mm"))

#Average fungi with functional annotation
p_FUNGI_func_ave_RA_BCNMDS <- plot_ordination(FUNGI_func_ave_RA, FUNGI_func_ave_RA_BCNMDS, type="samples", color="region", shape = "water_body") 
p_FUNGI_func_ave_RA_BCNMDS  <- p_FUNGI_func_ave_RA_BCNMDS + geom_point(size=8, alpha = 0.8) + ggtitle("Relative Abundance counts of Fungi-assigned ASVs averaged by water body", subtitle = "Bray-Curtis NMDS Stress = 0.220, ANOSIM R-Stat 0.7761, p value < 0.001")
p_FUNGI_func_ave_RA_BCNMDS $layers
p_FUNGI_func_ave_RA_BCNMDS $layers <- p_FUNGI_func_ave_RA_BCNMDS $layers[-1]
p_FUNGI_func_ave_RA_BCNMDS <- p_FUNGI_func_ave_RA_BCNMDS + scale_colour_manual(values = colour_palette)
p_FUNGI_func_ave_RA_BCNMDS + font_size

#use ggsave to save to your wd as a pdf and jpeg
ggsave("p_FUNGI_func_ave_RA_BCNMDS.jpg", width = 297, height = 210, units = c("mm"))
ggsave("p_FUNGI_func_ave_RA_BCNMDS.pdf", width = 297, height = 210, units = c("mm"))

#Average metazoa with functional annotation
p_METAZOA_func_ave_RA_BCNMDS <- plot_ordination(METAZOA_func_ave_RA, METAZOA_func_ave_RA_BCNMDS, type="samples", color="region", shape = "water_body") 
p_METAZOA_func_ave_RA_BCNMDS  <- p_METAZOA_func_ave_RA_BCNMDS+ geom_point(size=8, alpha = 0.8) + ggtitle("Relative Abundance counts of Metazoa-assigned ASVs averaged by water body", subtitle = "Bray-Curtis NMDS Stress = 0.167, ANOSIM R-Stat 0.649, p value < 0.001")
p_METAZOA_func_ave_RA_BCNMDS $layers
p_METAZOA_func_ave_RA_BCNMDS $layers <- p_METAZOA_func_ave_RA_BCNMDS $layers[-1]
p_METAZOA_func_ave_RA_BCNMDS <- p_METAZOA_func_ave_RA_BCNMDS + scale_colour_manual(values = colour_palette)
p_METAZOA_func_ave_RA_BCNMDS + font_size

#use ggsave to save to your wd as a pdf and jpeg
ggsave("p_METAZOA_func_ave_RA_BCNMDS.jpg", width = 297, height = 210, units = c("mm"))
ggsave("p_METAZOA_func_ave_RA_BCNMDS.pdf", width = 297, height = 210, units = c("mm"))

#Average protist with functional annotation
p_PROTIST_func_ave_RA_BCNMDS <- plot_ordination(PROTIST_func_ave_RA, PROTIST_func_ave_RA_BCNMDS, type="samples", color="region", shape = "water_body") 
p_PROTIST_func_ave_RA_BCNMDS  <- p_PROTIST_func_ave_RA_BCNMDS + geom_point(size=8, alpha = 0.8) + ggtitle("Relative Abundance counts of Protist-assigned ASVs averaged by water body", subtitle = "Bray-Curtis NMDS Stress = 0.192, ANOSIM R-Stat 0.6905, p value < 0.001")
p_PROTIST_func_ave_RA_BCNMDS $layers
p_PROTIST_func_ave_RA_BCNMDS $layers <- p_PROTIST_func_ave_RA_BCNMDS $layers[-1]
p_PROTIST_func_ave_RA_BCNMDS <- p_PROTIST_func_ave_RA_BCNMDS + scale_colour_manual(values = colour_palette) 
p_PROTIST_func_ave_RA_BCNMDS + font_size

#use ggsave to save to your wd as a pdf and jpeg
ggsave("p_PROTIST_func_ave_RA_BCNMDS.jpg", width = 297, height = 210, units = c("mm"))
ggsave("p_PROTIST_func_ave_RA_BCNMDS.pdf", width = 297, height = 210, units = c("mm"))

#Plotting your ordination - 16S ####
#Average 16S
p_ps_16S_ave_RA_BCNMDS <- plot_ordination(ps_16S_ave_RA, ps_16S_ave_RA_BCNMDS, type="samples", color="region", shape = "water_body") 
p_ps_16S_ave_RA_BCNMDS  <- p_ps_16S_ave_RA_BCNMDS + geom_point(size=8, alpha = 0.8) + ggtitle("Relative Abundance counts of 16S ASVs averaged by water body", subtitle = "Bray-Curtis NMDS Stress = 0.104, ANOSIM R-Stat 0.6725, p value < 0.001")
p_ps_16S_ave_RA_BCNMDS $layers
p_ps_16S_ave_RA_BCNMDS $layers <- p_ps_16S_ave_RA_BCNMDS $layers[-1]
p_ps_16S_ave_RA_BCNMDS <- p_ps_16S_ave_RA_BCNMDS + scale_colour_manual(values = colour_palette) 
p_ps_16S_ave_RA_BCNMDS + font_size

#use ggsave to save to your wd as a pdf and jpeg
ggsave("p_ps_16S_ave_RA_BCNMDS.jpg", width = 297, height = 210, units = c("mm"))
ggsave("p_ps_16S_ave_RA_BCNMDS.pdf", width = 297, height = 210, units = c("mm"))

#Average cyanobacteria
p_cyanobacteria_ave_RA_BCNMDS <- plot_ordination(cyanobacteria_ave_RA, cyanobacteria_ave_RA_BCNMDS, type="samples", color="region", shape = "water_body") 
p_cyanobacteria_ave_RA_BCNMDS  <- p_cyanobacteria_ave_RA_BCNMDS + geom_point(size=8, alpha = 0.8) + ggtitle("Relative Abundance counts of cyanobacteria-assigned ASVs averaged by water body", subtitle = "Bray-Curtis NMDS Stress = 0.171, ANOSIM R-Stat 0.5685, p value < 0.001")
p_cyanobacteria_ave_RA_BCNMDS $layers
p_cyanobacteria_ave_RA_BCNMDS$layers <- p_cyanobacteria_ave_RA_BCNMDS $layers[-1]
p_cyanobacteria_ave_RA_BCNMDS <- p_cyanobacteria_ave_RA_BCNMDS + scale_colour_manual(values = colour_palette) 
p_cyanobacteria_ave_RA_BCNMDS + font_size

#use ggsave to save to your wd as a pdf and jpeg
ggsave("p_cyanobacteria_ave_RA_BCNMDS.jpg", width = 297, height = 210, units = c("mm"))
ggsave("p_cyanobacteria_ave_RA_BCNMDS.pdf", width = 297, height = 210, units = c("mm"))

#Average proteobacteria
p_proteobacteria_ave_RA_BCNMDS <- plot_ordination(proteobacteria_ave_RA, proteobacteria_ave_RA_BCNMDS, type="samples", color="region", shape = "water_body") 
p_proteobacteria_ave_RA_BCNMDS  <- p_proteobacteria_ave_RA_BCNMDS + geom_point(size=8, alpha = 0.8) + ggtitle("Relative Abundance counts of proteobacteria-assigned ASVs averaged by water body", subtitle = "Bray-Curtis NMDS Stress = 0.105, ANOSIM R-Stat 0.636, p value < 0.001")
p_proteobacteria_ave_RA_BCNMDS $layers
p_proteobacteria_ave_RA_BCNMDS$layers <- p_proteobacteria_ave_RA_BCNMDS $layers[-1]
p_proteobacteria_ave_RA_BCNMDS <- p_proteobacteria_ave_RA_BCNMDS + scale_colour_manual(values = colour_palette) 
p_proteobacteria_ave_RA_BCNMDS + font_size

#use ggsave to save to your wd as a pdf and jpeg
ggsave("p_proteobacteria_ave_RA_BCNMDS.jpg", width = 297, height = 210, units = c("mm"))
ggsave("p_proteobacteria_ave_RA_BCNMDS.pdf", width = 297, height = 210, units = c("mm"))

#Anosim calculations####
#check for significant difference between regions with anosim

#performs a non-parametric test of the significance of the sample-grouping you provide against a permutation-based null distribution generated by randomly permuting the sample labels many times (often 999).

#Why ANOSIM?
#The ANOSIM test is similar to ANOVA hypothesis test, but it uses a dissimilarity matrix as input instead of raw data
#It is also non-parametric which means it does not assume much about our data, so it's a good bet for skewed microbial abundance data
#ANOSIM uses ranked dissimilarities instead of actual distances. 
#In our case, we want to see if there is a significant difference in the mocribla community composition of groups betweem sample regions
#If we see a difference, this suggests there is a distinction between areas

#ANOSIM - 18S ####
#Step 1. Prepare your data
#Vegan requires we load data with columns as OTUS and rows as samples
#We can pull a relative abundance table from the relative abundance transformed data above 

count_18S_ra <- abundances(ps_func_18S_ave_RA) 
count_18S_ra <- as.data.frame(count_18S_ra)
#This needs to be transposed so the samples are rows
count_18S_ra <- t(count_18S_ra)
rownames(count_18S_ra) # check your samples are rows

#Our NMDS aimed to identify differences between regions
#extract the region variable from the phyloseq object
region = get_variable(ps_func_18S_ave_RA, "region")

#Step 2. ANOSIM of 18S data
ano_18S = anosim(count_18S_ra, region, distance = "bray", permutations = 9999)
ano_18S
#ANOSIM statistic R: 0.7046
#Significance: 0.0001
#Permutation: free
#Number of permutations: 999

summary(ano_18S)
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.101 0.136 0.167 0.198 

#Dissimilarity ranks between and within classes:
                #  0%    25%   50%    75% 100%   N
#Between            23 167.50 283.5 390.25  496 408
#1_kuujjuarapik      1  23.25  71.0 108.25  298  28
#2_umiujaq          33 108.50 199.5 286.00  447  28
#3_cambridge_bay     3  13.50  25.5  42.25   96  28
#4_bylot_island     44  44.00  44.0  44.00   44   1
#5_resolute         62  62.00  62.0  62.00   62   1
#6_ellesmere_island  2   2.00   2.0   2.00    2   1
#7_ice_shelves      17  17.00  17.0  17.00   17   1 #this shows how the sample sites with two sample points don't have any difference between as there is only one 'in between'

#ANOSIM - fungi ####
#Step 1. Prepare your data
#Vegan requires we load data with columns as OTUS and rows as samples
#We can pull a relative abundance table from the relative abundance transformed data above 

count_fungi_ra <- abundances(FUNGI_func_ave_RA)
count_fungi_ra <- as.data.frame(count_fungi_ra)
#This needs to be transposed so the samples are rows
count_fungi_ra <- t(count_fungi_ra)
rownames(count_fungi_ra) # check your samples are rows

#Our NMDS aimed to identify differences between regions
#extract the region variable from the phyloseq object
region = get_variable(FUNGI_func_ave_RA, "region")

#Step 2. ANOSIM of fungi data
ano_fungi = anosim(count_fungi_ra, region, distance = "bray", permutations = 9999)
ano_fungi
#ANOSIM statistic R: 0.7761
#Significance: 0.0001
#Permutation: free
#Number of permutations: 9999

summary(ano_fungi)
#Upper quantiles of permutations (null model):
#90%    95%  97.5%    99% 
#0.0986 0.1313 0.1622 0.1984 

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0984 0.1305 0.1637 0.1944 

#Dissimilarity ranks between and within classes:
#  0%    25%   50%    75% 100%   N
#Between             11 181.00 291.5 394.25  481 408
#1_kuujjuarapik      21  71.50 117.5 155.50  299  28
#2_umiujaq           14  42.75  94.5 154.25  267  28
#3_cambridge_bay      1   7.75  27.0  60.50  164  28
#4_bylot_island      63  63.00  63.0  63.00   63   1
#5_resolute          92  92.00  92.0  92.00   92   1
#6_ellesmere_island 208 208.00 208.0 208.00  208   1
#7_ice_shelves       15  15.00  15.0  15.00   15   1

#ANOSIM - metazoa ####
#Step 1. Prepare your data
#Vegan requires we load data with columns as OTUS and rows as samples
#We can pull a relative abundance table from the relative abundance transformed data above 

count_metazoa_ra <- abundances(METAZOA_func_ave_RA)
count_metazoa_ra <- as.data.frame(count_metazoa_ra)
#This needs to be transposed so the samples are rows
count_metazoa_ra <- t(count_metazoa_ra)
rownames(count_metazoa_ra) # check your samples are rows

#Our NMDS aimed to identify differences between regions
#extract the region variable from the phyloseq object
region = get_variable(METAZOA_func_ave_RA, "region")

#Step 2. ANOSIM of metazoa data
ano_metazoa = anosim(count_metazoa_ra, region, distance = "bray", permutations = 9999)
ano_metazoa
#ANOSIM statistic R: 0.649
#Significance: 0.0001
#Permutation: free
#Number of permutations: 9999

summary(ano_metazoa)
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.106 0.140 0.167 0.202 

#Dissimilarity ranks between and within classes:
#  0%    25%   50%    75%  100%   N
#Between             15 166.75 282.5 387.25 491.5 408
#1_kuujjuarapik       3  12.75  33.5  79.75 190.0  28
#2_umiujaq           31 152.75 255.5 374.00 491.5  28
#3_cambridge_bay      1  17.75  27.5  46.25 162.0  28
#4_bylot_island      57  57.00  57.0  57.00  57.0   1
#5_resolute          62  62.00  62.0  62.00  62.0   1
#6_ellesmere_island 188 188.00 188.0 188.00 188.0   1
#7_ice_shelves       77  77.00  77.0  77.00  77.0   1

#ANOSIM - Protists ####
#Step 1. Prepare your data
#Vegan requires we load data with columns as OTUS and rows as samples
#We can pull a relative abundance table from the relative abundance transformed data above 

count_protist_ra <- abundances(PROTIST_func_ave_RA)
count_protist_ra <- as.data.frame(count_protist_ra)
#This needs to be transposed so the samples are rows
count_protist_ra <- t(count_protist_ra)
rownames(count_protist_ra) # check your samples are rows

#Our NMDS aimed to identify differences between regions
#extract the region variable from the phyloseq object
region = get_variable(PROTIST_func_ave_RA, "region")

#Step 2. ANOSIM of protist data
ano_protist = anosim(count_protist_ra, region, distance = "bray", permutations = 9999)
ano_protist
#ANOSIM statistic R: 0.6905
#Significance: 0.0001
#Permutation: free
#Number of permutations: 9999

summary(ano_protist)
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.102 0.135 0.162 0.190 

#Dissimilarity ranks between and within classes:
#  0%    25%   50%    75% 100%   N
#Between             9 171.50 284.5 390.50  496 408
#1_kuujjuarapik      1  51.00  95.5 156.75  418  28
#2_umiujaq          24  84.25 173.0 247.25  427  28
#3_cambridge_bay     3  16.75  32.0  44.75   97  28
#4_bylot_island     38  38.00  38.0  38.00   38   1
#5_resolute         86  86.00  86.0  86.00   86   1
#6_ellesmere_island  2   2.00   2.0   2.00    2   1
#7_ice_shelves      12  12.00  12.0  12.00   12   1

#ANOSIM - 16S ####
#Step 1. Prepare your data
#Vegan requires we load data with columns as OTUS and rows as samples
#We can pull a relative abundance table from the relative abundance transformed data above 

count_16S_ra <- abundances(ps_16S_ave_RA)
count_16S_ra <- as.data.frame(count_16S_ra)
#This needs to be transposed so the samples are rows
count_16S_ra <- t(count_16S_ra)
rownames(count_16S_ra) # check your samples are rows

#Our NMDS aimed to identify differences between regions
#extract the region variable from the phyloseq object
region = get_variable(ps_16S_ave_RA, "region")

#Step 2. ANOSIM of 16S data
ano_16S = anosim(count_16S_ra, region, distance = "bray", permutations = 9999)
ano_16S
#ANOSIM statistic R: 0.6725 
#Significance: 0.0001 
#Permutation: free
#Number of permutations: 9999

summary(ano_16S)
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.100 0.143 0.179 0.203 

#Dissimilarity ranks between and within Orderes:
#                    0%    25%   50%    75% 100%   N
#Between             33 171.75 282.5 388.25  495 408
#1_kuujjuarapik      11  41.75  74.0 142.00  320  28
#2_umiujaq            7 103.75 183.5 321.75  496  28
#3_cambridge_bay      1   8.75  17.5  28.25   61  28
#4_bylot_island      44  44.00  44.0  44.00   44   1
#5_resolute         122 122.00 122.0 122.00  122   1
#6_ellesmere_island  26  26.00  26.0  26.00   26   1
#7_ice_shelves       24  24.00  24.0  24.00   24   1 

#ANOSIM - cyanobacteria ####
#Step 1. Prepare your data
#Vegan requires we load data with columns as OTUS and rows as samples
#We can pull a relative abundance table from the relative abundance transformed data above 

count_cyanobacteria_ra <- abundances(cyanobacteria_ave_RA)
count_cyanobacteria_ra <- as.data.frame(count_cyanobacteria_ra)
#This needs to be transposed so the samples are rows
count_cyanobacteria_ra <- t(count_cyanobacteria_ra)
rownames(count_cyanobacteria_ra) # check your samples are rows

#Our NMDS aimed to identify differences between regions
#extract the region variable from the phyloseq object
region = get_variable(cyanobacteria_ave_RA, "region")

#Step 2. ANOSIM of cyanobacteria data
ano_cyanobacteria = anosim(count_cyanobacteria_ra, region, distance = "bray", permutations = 9999)
ano_cyanobacteria
#ANOSIM statistic R: 0.5685 
#Significance: 0.0001 

#Permutation: free
#Number of permutations: 9999

summary(ano_cyanobacteria)
#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0997 0.1316 0.1603 0.1966 

#Dissimilarity ranks between and within classes:
#  0%    25%   50%    75% 100%   N
#Between            28 164.75 275.5 385.25  495 408
#1_kuujjuarapik     18  49.50  95.5 122.25  312  28
#2_umiujaq          11 168.00 306.5 401.25  486  28
#3_cambridge_bay     1   7.75  16.5  38.00   76  28
#4_bylot_island     78  78.00  78.0  78.00   78   1
#5_resolute         77  77.00  77.0  77.00   77   1
#6_ellesmere_island 25  25.00  25.0  25.00   25   1
#7_ice_shelves      39  39.00  39.0  39.00   39   1

#ANOSIM - proteobacteria ####
#Step 1. Prepare your data
#Vegan requires we load data with columns as OTUS and rows as samples
#We can pull a relative abundance table from the relative abundance transformed data above 

count_proteobacteria_ra <- abundances(proteobacteria_ave_RA)
count_proteobacteria_ra <- as.data.frame(count_proteobacteria_ra)
#This needs to be transposed so the samples are rows
count_proteobacteria_ra <- t(count_proteobacteria_ra)
rownames(count_proteobacteria_ra) # check your samples are rows

#Our NMDS aimed to identify differences between regions
#extract the region variable from the phyloseq object
region = get_variable(proteobacteria_ave_RA,  "region")

#Step 2. ANOSIM of cyanobacteria data
ano_proteobacteria = anosim(count_proteobacteria_ra, region, distance = "bray", permutations = 9999)
ano_proteobacteria
#ANOSIM statistic R: 0.636 
#Significance: 0.0001 

#Permutation: free
#Number of permutations: 9999

summary(ano_proteobacteria)
#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0973 0.1305 0.1609 0.1998 

#Dissimilarity ranks between and within classes:
#  0%    25%   50%    75% 100%   N
#Between             23 169.75 281.5 386.25  496 408
#1_kuujjuarapik       9  45.00  81.5 140.00  302  28
#2_umiujaq            3 127.00 207.5 398.75  493  28
#3_cambridge_bay      1   9.50  19.5  30.25   70  28
#4_bylot_island      46  46.00  46.0  46.00   46   1
#5_resolute         119 119.00 119.0 119.00  119   1
#6_ellesmere_island  40  40.00  40.0  40.00   40   1
#7_ice_shelves       15  15.00  15.0  15.00   15   1

#CCA of Environmental Variables####
#1. Load your data
#16S count table
CCA_count_16S <- abundances(ps_16S_ave_RA)
CCA_count_16S <- as.data.frame(CCA_count_16S) #This needs to be transposed so the samples are rows for vegan
CCA_count_16S<- t(CCA_count_16S)
rownames(CCA_count_16S) #good are rows are now our samples

#16S metadata
CCA_meta_16S <- meta(ps_16S_ave_RA)
rownames(CCA_meta_16S)

#18S count table
CCA_count_18S <- abundances(ps_func_18S_ave_RA)
CCA_count_18S <- as.data.frame(CCA_count_18S) #This needs to be transposed so the samples are rows for vegan
CCA_count_18S<- t(CCA_count_18S)
rownames(CCA_count_18S) #good are rows are now our samples

#16S metadata
CCA_meta_18S <- meta(ps_func_18S_ave_RA)
rownames(CCA_meta_18S)

#2. 16S CCA####

#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
CCA_meta_16S<-CCA_meta_16S[rownames(CCA_count_16S),]
CCA_meta_16S

#Filter out any sample ASVs that have zero entries 
CCA_count_16S<-subset(CCA_count_16S,rowSums(CCA_count_16S)!=0)

#subset meta_table to desired environmental variables
colnames(CCA_meta_16S)
myvars <-c("water_temp", "water_pH","water_conductivity","ave_ann_temp","year")
myvars
CCA_meta_16S <- CCA_meta_16S[myvars]
colnames(CCA_meta_16S) 
#[1] "water_temp"         "water_pH"           "water_conductivity" "ave_ann_temp" "year" 

#we need listwise deletion of the temp samples that do not have data for shore pond 1-3
k <- complete.cases(CCA_meta_16S)
k

#adonis
#Info on https://chrischizinski.github.io/rstats/adonis/
#adonis works by first finding the centroids for each group and then calculates the squared deviations of each of site to that centroid. Then significance tests are performed using F-tests based on sequential sums of squares from permutations of the raw data.


#Use adonis to find significant environmental variables
CCA_count_16S.adonis <- adonis(CCA_count_16S[k,] ~ ., data=CCA_meta_16S[k,])
CCA_count_16S.adonis
#OUTPUT:
#                   Df  SumsOfSqs MeanSqs  F.Model R2       Pr(>F)    
# water_temp          1    0.8808 0.88075  3.3622 0.08042  0.001 ***
# water_pH            1    1.5077 1.50771  5.7556 0.13767  0.001 ***
# water_conductivity  1    0.5681 0.56814  2.1689 0.05188  0.004 ** 
# ave_ann_temp        1    0.5464 0.54642  2.0859 0.04989  0.016 *  
# sample_year_temp    1    0.6707 0.67071  2.5604 0.06124  0.006 ** 
# sample_month_temp   1    0.3441 0.34406  1.3134 0.03142  0.148    
# year                1    0.4092 0.40918  1.5620 0.03736  0.079 .  
# Residuals          23    6.0250 0.26195         0.55013           
# Total              30   10.9519                 1.00000           
# ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#INTERPRETATION: This shows us that water temperature and pH have a highly significant influence on ASV counts (0.00), conductivity and sample year temp have a significant impact (0.01) and there is a significant influence of average annual temperature (0.05)

#Extract the best variables (p <= 0.01)
bestEnvVariables<-rownames(CCA_count_16S.adonis$aov.tab)[CCA_count_16S.adonis$aov.tab$"Pr(>F)"<=0.01] #anything above 0.001
bestEnvVariables

#Last two are NA entries, so we have to remove them
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables

#We now use these significant values in our cca
eval(parse(text=paste("sol <- cca(CCA_count_16S[k,] ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=CCA_meta_16S[k,])",sep=""))) #this makes a sol variable

scrs<-scores(sol, display=c("sp","wa","lc","bp","cn"))
attributes(scrs)

#extract site data
df_sites<-data.frame(scrs$sites,t(as.data.frame(strsplit(rownames(scrs$sites),"_"))))
colnames(df_sites)<-c("x","y","SampleID")
df_sites

#add row with region info and pond info to your CCA dataframe
cbind_meta_16S <- meta(ps_16S_ave_RA)
cbind_meta_16S

#add row with region info and pond info to your CCA dataframe
pond_df = cbind_meta_16S$new_pond
pond_df = as.data.frame(pond_df, row.names = NULL)
pond_df = subset(pond_df, pond_df!="RE1") #RE1 doesn't have water data so had to remove it from CCA analysis
pond_df


region_df = cbind_meta_16S$region
region_df = as.data.frame(region_df, row.names=NULL)
region_df #we need to remove row 9 for the Resolute sample
region_df = region_df[-c(9), ]
region_df # should only be one resolute
water_body_df = cbind_meta_16S$water_body
water_body_df = as.data.frame(water_body_df, row.names = NULL)
water_body_df
water_body_df = water_body_df[-c(9), ]
df_sites <- cbind(df_sites,pond_df, region_df, water_body_df)
df_sites #we have now added our new sample labels and region metadata to our CCA plot

p<-ggplot()
p<-p+geom_point(data=df_sites, aes(x,y,colour=region_df, shape = water_body_df), alpha=0.8, size = 8)

#draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)

# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines 
# which are the coordinates of the heads of unit length vectors. In plot these are 
# scaled by their correlation (square root of the column r2) so that "weak" predictors 
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths 
# using command scores. The plotted (and scaled) arrows are further adjusted to the 
# current graph using a constant multiplier: this will keep the relative r2-scaled 
# lengths of the arrows but tries to fill the current plot. You can see the multiplier 
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul. 

#add df arrows
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)
df_arrows

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.9)
CCA_16S<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.8, size = 7) + scale_colour_manual(values = colour_palette) 
CCA_16S + font_size

#use ggsave to save to your wd as a pdf and jpeg
ggsave("CCA_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("CCA_16S.pdf", width = 297, height = 210, units = c("mm"))

#3. 18S CCA####
#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
CCA_meta_18S<-CCA_meta_18S[rownames(CCA_count_18S),]
CCA_meta_18S

#Filter out any sample ASVs that have zero entries 
CCA_count_18S<-subset(CCA_count_18S,rowSums(CCA_count_18S)!=0)

#subset meta_table to desired environmental variables
colnames(CCA_meta_18S)
myvars <-c("water_temp", "water_pH","water_conductivity","ave_ann_temp")
myvars
CCA_meta_18S <- CCA_meta_18S[myvars]
colnames(CCA_meta_18S) 
#[1] "water_temp"         "water_pH"           "water_conductivity" "ave_ann_temp" "year" 

#we need listwise deletion of the temp samples that do not have data for shore pond 1-3
k <- complete.cases(CCA_meta_18S)
k

#Use adonis to find significant environmental variables
CCA_count_18S.adonis <- adonis(CCA_count_18S[k,] ~ ., data=CCA_meta_18S[k,])
CCA_count_18S.adonis
#OUTPUT:
#                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#water_temp           1    0.9581 0.95805  2.8649 0.07739  0.001 ***
#water_pH             1    0.9663 0.96632  2.8896 0.07806  0.001 ***
#water_conductivity   1    0.6462 0.64621  1.9324 0.05220  0.001 ***
#ave_ann_temp         1    0.5394 0.53936  1.6129 0.04357  0.012 *  
#sample_year_temp     1    0.6561 0.65611  1.9620 0.05300  0.003 ** 
#sample_month_temp    1    0.4884 0.48841  1.4605 0.03945  0.036 *  
#year                 1    0.4337 0.43368  1.2969 0.03503  0.072 .  
#Residuals           23    7.6914 0.33441         0.62130           
#Total               30   12.3795                 1.00000           
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#INTERPRETATION: This shows us that water temperature, pH and conductivity have a highly significant influence on ASV counts (0.001), sample year temp has a significant impact (0.01) and there is a significant influence of average annual temperature (0.05)

#Extract the best variables (p <= 0.01)
bestEnvVariables<-rownames(CCA_count_18S.adonis$aov.tab)[CCA_count_18S.adonis$aov.tab$"Pr(>F)"<=0.01] #anything above 0.001
bestEnvVariables

#Last two are NA entries, so we have to remove them
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables

#We now use these significant values in our cca
eval(parse(text=paste("sol <- cca(CCA_count_18S[k,] ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=CCA_meta_18S[k,])",sep=""))) #this makes a sol variable

scrs<-scores(sol, display=c("sp","wa","lc","bp","cn"))
attributes(scrs)

#extract site data
df_sites<-data.frame(scrs$sites,t(as.data.frame(strsplit(rownames(scrs$sites),"_"))))
colnames(df_sites)<-c("x","y","SampleID")
df_sites

#add row with region info and pond info to your CCA dataframe
cbind_meta_18S <- meta(ps_func_18S_ave_RA)
cbind_meta_18S

#add row with region info and pond info to your CCA dataframe
pond_df = cbind_meta_18S$new_pond
pond_df = as.data.frame(pond_df, row.names = NULL)
pond_df = subset(pond_df, pond_df!="RE1") #RE1 doesn't have water data so had to remove it from CCA analysis
pond_df
region_df = cbind_meta_18S$region
region_df = as.data.frame(region_df, row.names=NULL)
region_df #we need to remove row 9 for the Resolute sample
region_df = region_df[-c(9), ]
region_df # should only be one resolute
water_body_df = cbind_meta_18S$water_body
water_body_df = as.data.frame(water_body_df, row.names = NULL)
water_body_df
water_body_df = water_body_df[-c(9), ]
df_sites <- cbind(df_sites,pond_df, region_df, water_body_df)
df_sites #we have now added our new sample labels and region metadata to our CCA plot

p<-ggplot()
p<-p+geom_point(data=df_sites, aes(x,y,colour=region_df, shape = water_body_df), alpha=0.8, size = 8)

#draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)

# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines 
# which are the coordinates of the heads of unit length vectors. In plot these are 
# scaled by their correlation (square root of the column r2) so that "weak" predictors 
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths 
# using command scores. The plotted (and scaled) arrows are further adjusted to the 
# current graph using a constant multiplier: this will keep the relative r2-scaled 
# lengths of the arrows but tries to fill the current plot. You can see the multiplier 
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul. 

#add df arrows
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)
df_arrows

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="#808080",alpha=0.9)
CCA_18S<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="#808080",alpha=0.8, size = 7) + scale_colour_manual(values = colour_palette)
CCA_18S + font_size

#use ggsave to save to your wd as a pdf and jpeg
ggsave("CCA_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("CCA_18S.pdf", width = 297, height = 210, units = c("mm"))

#adds ASVs to plot
#df_species<- as.data.frame(scrs$species)
#colnames(df_species)<-c("x","y")
#p<-p+geom_point(data=df_species,aes(x,y,shape="Species"), alpha=0.01)

#16S - INDICATOR SPECIES ANALYSIS - REGION ####
#When analysing microbial data you may want to identify the microbial species that are found more often in more treatment group than another. 
#There are a number of ways to look at this using statistical tests
#Indicator species in OTU tables are those that correlate with particular sample variables, in our case the regions

#There is a specific package for Indicator Species Analysis
#install.packages("indicspecies")
library(indicspecies)

#Agglomerate your phyloseq objects for 18S by order, NArm = FALSE (do not remove ASVs unassigned at order level)
ps_16S_ave_RA
unique(tax_table(ps_16S_ave_RA)[,"Order"] ) #this tells you that there are 292 unique orders
ps_16S_ave_RA_Order <- tax_glom(ps_16S_ave_RA, "Order", NArm = TRUE)
ps_16S_ave_RA_Order #this shows us all our ASVs have been sorted into their respective 381 taxa now, this is because tax glom will not merge at Order level if they have different higher level classification
unique(tax_table(ps_16S_ave_RA_Order)[,"Order"] )

#Now we replace the ASVIDs with the order, sadly we can't use this function at Order level as there are multiple orders with the same name
#taxa_names(ps_16S_ave_RA_Order) <- tax_table#(ps_16S_ave_RA_order)[,"Order"]
#taxa_names(ps_16S_ave_RA_Order)

#An alternative to the above using Microbiome utilities package
ps_16S_ave_RA_Order <- microbiomeutilities::format_to_besthit(ps_16S_ave_RA_Order)
head(otu_table(ps_16S_ave_RA_Order)) #you can now see they have the species names

#Make a count table
ps_16S_ave_RA_Order_count <- abundances(ps_16S_ave_RA_Order) #this utilises the compositional function in microbiome package to create relative abundances of our data
ps_16S_ave_RA_Order_count <- as.data.frame(ps_16S_ave_RA_Order_count)
#This needs to be transposed so the samples are rows for vegan
ps_16S_ave_RA_Order_count <- t(ps_16S_ave_RA_Order_count)
rownames(ps_16S_ave_RA_Order_count) #good are rows are now our samples

#We need to create a table with our region details from our metadata table
indic_region <- as.data.frame(meta_16S)
colnames(indic_region) #regions is column 3, water_body is column 5
indic_region = indic_region[,3:5]
indic_region 

#Calculate inv variable for 16S phyla
inv_region_16s <- multipatt(ps_16S_ave_RA_Order_count, indic_region$region, func = "r.g", control = how(nperm=9999))
summary(inv_region_16s)

#INDICATOR SPECIES RESULTS - REGION 16S####
#Total number of species: 291
#Selected number of species: 27 
#Number of species associated to 1 group: 17 
#Number of species associated to 2 groups: 5 
#Number of species associated to 3 groups: 2 
#Number of species associated to 4 groups: 3 
#Number of species associated to 5 groups: 0 
#Number of species associated to 6 groups: 0 

#List of species associated to each combination: 
  
#  Group 1_kuujjuarapik  #sps.  1 
#stat p.value  
#ASV3969:o__Candidatus_Moranbacteria 0.743  0.0276 *
  
#  Group 4_bylot_island  #sps.  10 
#stat p.value   
#ASV111:o__Chloroflexales                  0.839  0.0056 **
#  ASV164:o__Thermosynechococcales           0.835  0.0118 * 
#  ASV270:o__Defluviicoccales                0.800  0.0146 * 
#  ASV1467:o__Micromonosporales              0.759  0.0264 * 
#  ASV2355:o__Chthonomonadales               0.759  0.0246 * 
#  ASV1322:o__Micropepsales                  0.757  0.0235 * 
#  ASV2621:o__Elev-16S-1166                  0.752  0.0217 * 
#  ASV209:o__Pseudonocardiales               0.744  0.0170 * 
#  ASV536:o__Oxyphotobacteria_Incertae_Sedis 0.698  0.0400 * 
#  ASV1161:o__Thermomicrobiales              0.694  0.0359 * 
  
#  Group 5_resolute  #sps.  4 
#stat p.value  
#ASV2071:o__Desulfovibrionales 0.814  0.0230 *
#  ASV1018:o__mle1-8             0.791  0.0160 *
#  ASV3979:o__Vampirovibrionales 0.711  0.0322 *
#  ASV1609:o__Competibacterales  0.675  0.0345 *
#  
#  Group 6_ellesmere_island  #sps.  2 
#stat p.value  
#ASV4720:o__Zavarziniales                     0.724  0.0326 *
#  ASV511:o__Gammaproteobacteria_Incertae_Sedis 0.687  0.0440 *
#  
#  Group 3_cambridge_bay+7_ice_shelves  #sps.  1 
#stat p.value   
#ASV17:o__Phormidesmiales 0.814  0.0067 **
#  
#  Group 4_bylot_island+5_resolute  #sps.  1 
#stat p.value  
#ASV139:o__Haliangiales 0.741  0.0288 *
#  
# Group 4_bylot_island+6_ellesmere_island  #sps.  2 
#stat p.value  
#ASV5233:o__0319-6G20       0.699  0.0342 *
#ASV60:o__Corynebacteriales 0.688  0.0355 *
  
#  Group 5_resolute+7_ice_shelves  #sps.  1 
#stat p.value  
#ASV295:o__RD017 0.792  0.0206 *
  
#  Group 3_cambridge_bay+5_resolute+7_ice_shelves  #sps.  1 
#stat p.value  
#ASV100:o__SM1A07 0.786  0.0134 *
  
#  Group 4_bylot_island+5_resolute+6_ellesmere_island  #sps.  1 
#stat p.value  
#ASV208:o__Tistrellales 0.784  0.0107 *
  
#  Group 3_cambridge_bay+4_bylot_island+5_resolute+7_ice_shelves  #sps.  1 
#stat p.value  
#ASV95:o__Nannocystales 0.682  0.0334 *
  
#  Group 3_cambridge_bay+5_resolute+6_ellesmere_island+7_ice_shelves  #sps.  1 
#stat p.value    
#ASV303:o__Phycisphaerales 0.83   3e-04 ***
  
#  Group 4_bylot_island+5_resolute+6_ellesmere_island+7_ice_shelves  #sps.  1 
#stat p.value  
#ASV219:o__Gemmatales 0.699  0.0322 *
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#16S - INDICATOR SPECIES ANALYSIS - WATER BODY####
#When analysing microbial data you may want to identify the microbial species that are found more often in more treatment group than another. 
#There are a number of ways to look at this using statistical tests
#Indicator species in OTU tables are those that correlate with particular sample variables, this time we will look at the water bodies in our bray-curtis to see if we identify any distinction between water body type, especially with the distinction seen in the rock_pool


#Calculate inv variable for 16S orders by water body
inv_water_body_16S <- multipatt(ps_16S_ave_RA_Order_count, indic_region$water_body, func = "r.g", control = how(nperm=9999))
summary(inv_water_body_16S)

#INDICATOR SPECIES WATER BODY - 16S RESULTS####
Association function: r.g
Significance level (alpha): 0.05

Total number of species: 291
Selected number of species: 57 
Number of species associated to 1 group: 57 
Number of species associated to 2 groups: 0 
Number of species associated to 3 groups: 0 

List of species associated to each combination: 
  
# Group rock_pool  #sps.  57 
# stat p.value   
#ASV9749:o__Candidatus_Peribacteria      1.000  0.0293 * 
#  ASV10984:o__Candidatus_Colwellbacteria  1.000  0.0293 * 
#  ASV15123:o__Clostridia_vadinBB60_group  1.000  0.0293 * 
#  ASV6719:o__Clostridia_UCG-014           1.000  0.0158 * 
#  ASV1049:o__Cloacimonadales              1.000  0.0126 * 
#  ASV3271:o__Oscillospirales              0.997  0.0021 **
#  ASV1050:o__Aminicenantales              0.995  0.0117 * 
#  ASV6793:o__Victivallales                0.989  0.0061 **
#  ASV4693:o__Woesearchaeales              0.988  0.0122 * 
#  ASV1874:o__Syntrophobacterales          0.986  0.0044 **
#  ASV16953:o__FCPU453                     0.986  0.0216 * 
#  ASV3581:o__MSBL9                        0.985  0.0015 **
#  ASV15732:o__Micrarchaeales              0.983  0.0243 * 
#  ASV10718:o__Lineage_IV                  0.982  0.0112 * 
#  ASV171:o__Methanobacteriales            0.982  0.0069 **
#  ASV11827:o__Candidatus_Lloydbacteria    0.982  0.0044 **
#  ASV790:o__Acidobacteriales              0.979  0.0022 **
#  ASV8474:o__Candidatus_Nomurabacteria    0.977  0.0088 **
#  ASV567:o__Pedosphaerales                0.975  0.0034 **
#  ASV1981:o__Saccharimonadales            0.973  0.0029 **
#  ASV5966:o__Candidatus_Magasanikbacteria 0.971  0.0282 * 
#  ASV12559:o__DG-20                       0.963  0.0442 * 
#  ASV15252:o__MVP-88                      0.960  0.0134 * 
#  ASV8820:o__Defferrisomatales            0.960  0.0125 * 
#  ASV5410:o__Omnitrophales                0.958  0.0246 * 
#  ASV642:o__Methylococcales               0.956  0.0076 **
#  ASV16955:o__MVP-21                      0.955  0.0283 * 
#  ASV6483:o__WCHB1-41                     0.953  0.0339 * 
#  ASV12722:o__Endomicrobiales             0.952  0.0331 * 
#  ASV2303:o__Methanosarciniales           0.952  0.0188 * 
#  ASV5445:o__S-BQ2-57_soil_group          0.950  0.0232 * 
#  ASV98:o__Bacteroidales                  0.950  0.0284 * 
#  ASV162:o__Holophagales                  0.936  0.0282 * 
#  ASV7392:o__Brevinematales               0.931  0.0366 * 
#  ASV879:o__Kryptoniales                  0.925  0.0279 * 
#  ASV6080:o__Candidatus_Yanofskybacteria  0.924  0.0229 * 
#  ASV16366:o__Candidatus_Buchananbacteria 0.920  0.0405 * 
#  ASV12565:o__Candidatus_Shapirobacteria  0.913  0.0386 * 
#  ASV4976:o__Caldisericales               0.912  0.0320 * 
#  ASV3234:o__Kiritimatiellales            0.904  0.0295 * 
#  ASV77:o__Frankiales                     0.892  0.0080 **
#  ASV197:o__OPB41                         0.892  0.0207 * 
#  ASV5666:o__vadinBA26                    0.888  0.0346 * 
#  ASV2881:o__SJA-15                       0.883  0.0440 * 
#  ASV374:o__Isosphaerales                 0.877  0.0172 * 
#  ASV563:o__Anaerolineales                0.876  0.0304 * 
#  ASV6414:o__Oligosphaerales              0.875  0.0368 * 
#  ASV78:o__Myxococcales                   0.872  0.0027 **
#  ASV537:o__Spirochaetales                0.848  0.0318 * 
#  ASV8823:o__AKIW659                      0.837  0.0410 * 
#  ASV7314:o__Syntrophomonadales           0.829  0.0472 * 
#  ASV240:o__Sphingobacteriales            0.825  0.0422 * 
#  ASV10523:o__FW22                        0.817  0.0419 * 
#  ASV6473:o__Syntrophorhabdales           0.802  0.0327 * 
#  ASV1782:o__Leptospirales                0.772  0.0342 * 
#  ASV3379:o__Methylomirabilales           0.764  0.0311 * 
#  ASV244:o__Geobacterales                 0.729  0.0421 * 
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#INDICATOR SPECIES - REGION 18S####

#Agglomerate your phyloseq objects for 18S by order, NArm = FALSE (do not remove ASVs unassigned at order level)
ps_func_18S_ave_RA
unique(tax_table(ps_func_18S_ave_RA)[,"Order"] ) #this tells you that there are 214 unique orders including NAs
ps_func_18S_ave_RA_order <- tax_glom(ps_func_18S_ave_RA, "Order", NArm = TRUE)
ps_func_18S_ave_RA_order #this shows us all our ASVs have been sorted into their respective 271 taxa now, this is because tax glom will not merge at Order level if they have different higher level Orderification
unique(tax_table(ps_func_18S_ave_RA_order)[,"Order"] ) #213 unique orders now

#Now we replace the ASVIDs with the order, sadly we can't use this function at Order level as there are multiple orders with the same name
#taxa_names(ps_func_18S_ave_RA) <- tax_table(ps_func_18S_ave_RA)[,"Order"]
#taxa_names(ps_func_18S_ave_RA)

#An alternative to the above using Microbiome utilities package
tax_table(ps_func_18S_ave_RA_order) <- tax_table(ps_func_18S_ave_RA_order)[,1:7]
ps_func_18S_ave_RA_order <- microbiomeutilities::format_to_besthit(ps_func_18S_ave_RA_order)
head(otu_table(ps_func_18S_ave_RA_order)) #you can now see they have the species names

#Make a count table
ps_func_18S_ave_RA_order_count <- abundances(ps_func_18S_ave_RA_order) 
ps_func_18S_ave_RA_order_count <- as.data.frame(ps_func_18S_ave_RA_order_count)
#This needs to be transposed so the samples are rows for vegan
ps_func_18S_ave_RA_order_count <- t(ps_func_18S_ave_RA_order_count)
rownames(ps_func_18S_ave_RA_order_count) #good our rows are now our samples

#We need to create a table with our region details from our metadata table
indic_region_18S <- as.data.frame(meta_18S)
colnames(indic_region_18S) #regions is column 3, water_body is column 5
indic_region_18S = indic_region_18S[,3:5]
indic_region_18S

#Calculate inv variable for 18S orders by region
inv_region_18S <- multipatt(ps_func_18S_ave_RA_order_count, indic_region_18S$region, func = "r.g", control = how(nperm=9999))
summary(inv_region_18S)

#INDICATOR SPECIES REGION RESULTS 18S####
#Association function: r.g
#Significance level (alpha): 0.05

#Total number of species: 213
#Selected number of species: 17 
#Number of species associated to 1 group: 9 
#Number of species associated to 2 groups: 3 
#Number of species associated to 3 groups: 4 
#Number of species associated to 4 groups: 1 
#Number of species associated to 5 groups: 0 
#Number of species associated to 6 groups: 0 

#List of species associated to each combination: 
  
#  Group 4_bylot_island  #sps.  5 
#stat p.value   
#ASV1:f__Sphaeropleales      0.952  0.0028 **
#  ASV476:f__Hyphochytrydiales 0.892  0.0065 **
#  ASV414:f__Cryptofilida      0.813  0.0189 * 
#  ASV626:f__Lobosa_XX         0.772  0.0181 * 
#  ASV544:f__Variosea_X        0.737  0.0125 * 
#  
#  Group 5_resolute  #sps.  2 
#stat p.value   
#ASV3279:f__Gonyaulacales 0.783  0.0079 **
#  ASV6:f__Annelida_X       0.747  0.0348 * 
#  
#  Group 6_ellesmere_island  #sps.  2 
#stat p.value  
#ASV5131:f__Ustilaginomycotina 0.801  0.0210 *
#  ASV68:f__Saccharomycotina     0.693  0.0427 *
#  
#  Group 4_bylot_island+5_resolute  #sps.  1 
#stat p.value  
#ASV2591:f__Haptoria_5 0.691  0.0279 *
  
#  Group 5_resolute+6_ellesmere_island  #sps.  2 
#stat p.value  
#ASV2430:f__Phaeophyceae_X 0.741  0.0252 *
#  ASV1131:f__Plagiopylea_X  0.710  0.0425 *
#  
#  Group 3_cambridge_bay+5_resolute+7_ice_shelves  #sps.  1 
#stat p.value  
#ASV10:f__Turbellaria 0.735  0.0274 *
#  
#  Group 4_bylot_island+5_resolute+6_ellesmere_island  #sps.  1 
#stat p.value  
#ASV721:f__Prokinetoplastida 0.734  0.0455 *
#  
#  Group 5_resolute+6_ellesmere_island+7_ice_shelves  #sps.  2 
#stat p.value  
#ASV272:f__Gymnodiniales 0.775  0.0118 *
#  ASV52:f__Gastrotricha_X 0.723  0.0320 *
#  
#  Group 3_cambridge_bay+5_resolute+6_ellesmere_island+7_ice_shelves  #sps.  1 
#stat p.value   
#ASV3:f__Chromadorea 0.764   0.003 **
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#

##18S - INDICATOR SPECIES ANALYSIS - WATER BODY####
#Calculate inv variable for 18S orders by water body
inv_water_body_18S <- multipatt(ps_func_18S_ave_RA_order_count, indic_region_18S$water_body, func = "r.g", control = how(nperm=9999))
summary(inv_water_body_18S)

#INDICATOR SPECIES RESULTS WATER BODY 18S####
#Total number of species: 213
#Selected number of species: 8 
#Number of species associated to 1 group: 7 
#Number of species associated to 2 groups: 1 
#Number of species associated to 3 groups: 0 

#List of species associated to each combination: 
  
#  Group rock_pool  #sps.  7 
#                                 stat p.value   
#  ASV5120:f__Peritrichia_1       0.978  0.0188 * 
#  ASV2182:f__Tetramitia_III      0.932  0.0198 * 
#  ASV1595:f__Tectofilosida       0.931  0.0115 * 
#  ASV40:f__Eustigmatophyceae_X   0.925  0.0052 **
#  ASV3600:f__Monomastigales      0.907  0.0012 **
#  ASV3573:f__Oligohymenophorea_X 0.885  0.0423 * 
#  ASV2754:f__Heterolobosea_X     0.753  0.0478 * 
  
#  Group lake+meltwater_pond  #sps.  1 
#                     stat    p.value  
#ASV3:f__Chromadorea 0.747   0.047 *

#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

##16S - SIMPER ANALYSIS####
#Simper produces a lot of results as it compares all your ASVS, if we did that with our dataset, we'd have 15250 comparisons, too many!
#Let's use phyla for simper analysis, as we want to see the contribution of the different phyla to the regions
#Agglomerate your phyloseq objects for 16S by phylum, NArm = TRUE (we will remove ASVs unassigned at phylum level)
ps_16S_ave_RA
unique(tax_table(ps_16S_ave_RA)[,"Order"] ) #this tells you that there are 49 unique phyla
ps_16S_ave_RA_Order <- tax_glom(ps_16S_ave_RA, "Order", NArm = TRUE)
ps_16S_ave_RA_Order #this shows us all our ASVs have been sorted into their respective 49 taxa now, this is because tax glom will not merge at Order level if they have different higher level Orderification
unique(tax_table(ps_16S_ave_RA_Phylum)[,"Phylum"] )

#Now we replace the ASVIDs with the Division names!
taxa_names(ps_16S_ave_RA_Phylum) <- tax_table(ps_16S_ave_RA_Phylum)[,"Phylum"]
taxa_names(ps_16S_ave_RA_Phylum)

#Make a count table
ps_16S_ave_RA_Phylum_count <- abundances(ps_16S_ave_RA_Phylum) #this utilises the compositional function in microbiome package to create relative abundances of our data
ps_16S_ave_RA_Phylum_count <- as.data.frame(ps_16S_ave_RA_Phylum_count)
#This needs to be transposed so the samples are rows for vegan
ps_16S_ave_RA_Phylum_count <- t(ps_16S_ave_RA_Phylum_count)
rownames(ps_16S_ave_RA_Phylum_count) #good are rows are now our samples

#We need to create a table with our region details from our metadata table
simper_meta_16S <- as.data.frame(meta_16S)
colnames(simper_meta_16S) #regions is column 3, water_body is column 5
simper_meta_16S = simper_meta_16S[,3:5]
simper_meta_16S

#SIMPER 16S REGION####
simp_16S_region <- simper(ps_16S_ave_RA_Order_count, simper_meta_16S$region, permutations = 9999)
simp_16S_region

#SIMPER 16S REGION RESULTS####
#cumulative contributions of most influential species:
  
#  $`1_kuujjuarapik_2_umiujaq`
#Proteobacteria  Actinobacteriota     Cyanobacteria      Bacteroidota Verrucomicrobiota       Chloroflexi 
#0.1776131         0.3418348         0.4818256         0.5908437         0.6649919         0.7389898 

#$`1_kuujjuarapik_3_cambridge_bay`
#Actinobacteriota    Cyanobacteria   Proteobacteria     Bacteroidota      Chloroflexi 
#0.2156193        0.3981251        0.5458090        0.6437111        0.7037836 

#$`1_kuujjuarapik_4_bylot_island`
#Actinobacteriota   Proteobacteria      Chloroflexi    Cyanobacteria  Planctomycetota      Myxococcota 
#0.1660754        0.3211825        0.4659507        0.6029951        0.6979063        0.7685284 

#$`1_kuujjuarapik_5_resolute`
#Actinobacteriota    Cyanobacteria   Proteobacteria  Planctomycetota     Bacteroidota 
#0.2345092        0.4258853        0.5883795        0.6751119        0.7434938 

#$`1_kuujjuarapik_6_ellesmere_island`
#Actinobacteriota    Cyanobacteria   Proteobacteria  Planctomycetota       Firmicutes 
#0.1887110        0.3627608        0.4956585        0.6281727        0.7162943 

#$`1_kuujjuarapik_7_ice_shelves`
#Actinobacteriota    Cyanobacteria   Proteobacteria     Bacteroidota  Planctomycetota 
#0.2071131        0.4056610        0.5743883        0.6971132        0.7847865 

#$`2_umiujaq_3_cambridge_bay`
#Cyanobacteria    Proteobacteria      Bacteroidota  Actinobacteriota       Chloroflexi Verrucomicrobiota   Acidobacteriota 
#0.1805344         0.3459755         0.4513130         0.5474575         0.6255920         0.6982311         0.7490499 

#$`2_umiujaq_4_bylot_island`
#Proteobacteria     Cyanobacteria      Bacteroidota       Chloroflexi  Actinobacteriota Verrucomicrobiota   Planctomycetota 
#0.1507734         0.2859411         0.3971149         0.5041755         0.5910345         0.6750363         0.7514897 

#$`2_umiujaq_5_resolute`
#Cyanobacteria    Proteobacteria  Actinobacteriota      Bacteroidota Verrucomicrobiota   Planctomycetota       Chloroflexi 
#0.2033570         0.3500471         0.4621943         0.5569725         0.6290624         0.6972923         0.7575557 

#$`2_umiujaq_6_ellesmere_island`
#Cyanobacteria   Proteobacteria     Bacteroidota  Planctomycetota       Firmicutes Actinobacteriota      Chloroflexi 
#0.1478526        0.2901538        0.4219337        0.5199883        0.6026656        0.6752494        0.7446765 

#$`2_umiujaq_7_ice_shelves`
#Cyanobacteria    Proteobacteria      Bacteroidota  Actinobacteriota   Planctomycetota Verrucomicrobiota 
#0.2215135         0.3633247         0.4798033         0.5803796         0.6569005         0.7267702 

#$`3_cambridge_bay_4_bylot_island`
#Cyanobacteria       Chloroflexi    Proteobacteria  Actinobacteriota      Bacteroidota Verrucomicrobiota 
#0.1555505         0.2975641         0.4314187         0.5583723         0.6654611         0.7456708 

#$`3_cambridge_bay_5_resolute`
#Proteobacteria    Cyanobacteria     Bacteroidota      Chloroflexi       Firmicutes  Planctomycetota Actinobacteriota      Myxococcota 
#0.1888416        0.3199143        0.4215518        0.5172726        0.5788408        0.6387525        0.6972382        0.7529355 

#$`3_cambridge_bay_6_ellesmere_island`
#Bacteroidota   Cyanobacteria  Proteobacteria      Firmicutes Planctomycetota     Chloroflexi 
#0.1718797       0.3205491       0.4536951       0.5739274       0.6636166       0.7330439 

#$`3_cambridge_bay_7_ice_shelves`
#Cyanobacteria  Proteobacteria    Bacteroidota     Chloroflexi Planctomycetota 
#0.2097406       0.3891084       0.5571430       0.6366361       0.7105073 

#$`4_bylot_island_5_resolute`
#Cyanobacteria  Actinobacteriota       Chloroflexi Verrucomicrobiota 
#0.2104808         0.4074988         0.5956611         0.7125518 

#`4_bylot_island_6_ellesmere_island`
#Chloroflexi     Cyanobacteria        Firmicutes Verrucomicrobiota  Actinobacteriota       Myxococcota 
#0.1851996         0.3297660         0.4440200         0.5422707         0.6359065         0.7121751 

#$`4_bylot_island_7_ice_shelves`
#Cyanobacteria      Bacteroidota       Chloroflexi  Actinobacteriota Verrucomicrobiota 
#0.1998404         0.3549777         0.4993228         0.6333682         0.7247421 

#$`5_resolute_6_ellesmere_island`
#Bacteroidota    Cyanobacteria       Firmicutes Actinobacteriota   Proteobacteria  Planctomycetota 
#0.1645438        0.3243969        0.4672482        0.5827136        0.6871839        0.7603129 

#$`5_resolute_7_ice_shelves`
#Cyanobacteria    Bacteroidota  Proteobacteria Planctomycetota     Myxococcota Gemmatimonadota 
#0.2374274       0.4450521       0.5255842       0.5942641       0.6570553       0.7046636 

#$`6_ellesmere_island_7_ice_shelves`
#Bacteroidota    Cyanobacteria   Proteobacteria       Firmicutes Actinobacteriota 
#0.2297084        0.4280769        0.5545486        0.6594165        0.7450808 


#SIMPER 16S WATER BODY####
simp_16S_water_body <- simper(ps_16S_ave_RA_Phylum_count, simper_meta_16S$water_body, permutations = 9999)
simp_16S_water_body

#SIMPER 16S WATER BODY RESULTS####
#cumulative contributions of most influential species:

#$pond_lake
#Cyanobacteria   Proteobacteria Actinobacteriota     Bacteroidota      Chloroflexi  Planctomycetota 
#0.1808445        0.3430216        0.4729231        0.5881595        0.6674008        0.7355569 

#$pond_rock_pool
#Cyanobacteria  Actinobacteriota    Proteobacteria      Bacteroidota   Acidobacteriota Verrucomicrobiota       Chloroflexi 
#0.1965707         0.3229897         0.4432563         0.5175301         0.5826198         0.6414909         0.6958778 
#Planctomycetota 
#0.7397798 

#pond_meltwater_pond
#Cyanobacteria    Proteobacteria      Bacteroidota  Actinobacteriota   Planctomycetota Verrucomicrobiota 
#0.2069180         0.3502832         0.4882084         0.6221982         0.6946932         0.7567242 

#$lake_rock_pool
#Cyanobacteria Actinobacteriota   Proteobacteria  Acidobacteriota     Bacteroidota  Planctomycetota      Chloroflexi 
#0.2694241        0.3917757        0.5034599        0.5787317        0.6340419        0.6824903        0.7306453 

#$lake_meltwater_pond
#Proteobacteria    Cyanobacteria     Bacteroidota  Planctomycetota Actinobacteriota 
#0.2209867        0.4149382        0.5859375        0.6827339        0.7510366 

#$rock_pool_meltwater_pond
#Cyanobacteria Actinobacteriota     Bacteroidota  Acidobacteriota  Planctomycetota 
#0.3071103        0.4643668        0.5514831        0.6293222        0.7058986 



#18S SIMPER ANALYSIS ####

#Agglomerate your phyloseq objects for 18S by Division, NArm = FALSE (do not remove ASVs unassigned at order level)
ps_func_18S_ave_RA
unique(tax_table(ps_func_18S_ave_RA)[,"Division"] ) #this tells you that there are 30 unique divisions
ps_18S_ave_RA_Division <- tax_glom(ps_func_18S_ave_RA, "Division", NArm = TRUE) # we choose true here to get rid of the NA assigned so we can focus on the interesting groups 
ps_18S_ave_RA_Division #this shows us all our ASVs have been sorted into 29 taxa, this is because we deterimined that NA was true this time so there is no NA
unique(tax_table(ps_18S_ave_RA_Division)[,"Division"] )

#Now we replace the ASVIDs with the Division names, we can't do this due to shared 
taxa_names(ps_18S_ave_RA_Division) <- tax_table(ps_18S_ave_RA_Division)[,"Division"]
taxa_names(ps_18S_ave_RA_Division)

#Make a count table
ps_18S_ave_RA_Division_count <- abundances(ps_18S_ave_RA_Division) #this utilises the compositional function in microbiome package to create relative abundances of our data
ps_18S_ave_RA_Division_count <- as.data.frame(ps_18S_ave_RA_Division_count)
#This needs to be transposed so the samples are rows for vegan
ps_18S_ave_RA_Division_count <- t(ps_18S_ave_RA_Division_count)
rownames(ps_18S_ave_RA_Division_count) #good are rows are now our samples

#We need to create a table with our region details from our metadata table
simper_meta_18S <- as.data.frame(meta_18S)
colnames(simper_meta_18s) #regions is column 3, water_body is column 5
simper_meta_18S  = simper_meta_18S [,3:5]
simper_meta_18S  

#region simper####
simp_18S_region <- simper(ps_18S_ave_RA_Division_count, simper_meta_18S$region, permutations = 9999)
simp_18S_region

#SIMPER 18s REGION OUTPUT####
#cumulative contributions of most influential species:
  
#  $`1_kuujjuarapik_2_umiujaq`
#Metazoa Dinoflagellata    Chlorophyta     Ochrophyta 
#0.2425791      0.4220428      0.5655863      0.7073431 

#$`1_kuujjuarapik_3_cambridge_bay`
#Metazoa    Chlorophyta     Ochrophyta Dinoflagellata 
#0.3296459      0.4844441      0.6328802      0.7773977 

#$`1_kuujjuarapik_4_bylot_island`
#Chlorophyta        Metazoa Dinoflagellata 
#0.4091666      0.5776690      0.7074414 

#$`1_kuujjuarapik_5_resolute`
#Metazoa    Chlorophyta Dinoflagellata     Ochrophyta 
#0.4035025      0.5545443      0.6938936      0.8056090 

#$`1_kuujjuarapik_6_ellesmere_island`
#Metazoa     Ochrophyta Dinoflagellata    Chlorophyta 
#0.2665193      0.4553872      0.6075914      0.7419297 

#$`1_kuujjuarapik_7_ice_shelves`
#Metazoa     Ochrophyta    Chlorophyta Dinoflagellata 
#0.3417118      0.5091095      0.6674780      0.8148899 

#$`2_umiujaq_3_cambridge_bay`
#Metazoa     Ochrophyta    Chlorophyta Dinoflagellata          Fungi 
#0.2932076      0.4445127      0.5755434      0.6669887      0.7384325 

#$`2_umiujaq_4_bylot_island`
#Chlorophyta     Metazoa  Ochrophyta 
#0.4495824   0.6246515   0.7151126 

#$`2_umiujaq_5_resolute`
#Metazoa    Chlorophyta     Ochrophyta Dinoflagellata 
#0.3875374      0.5173428      0.6359028      0.7228672 

#$`2_umiujaq_6_ellesmere_island`
#Metazoa     Ochrophyta          Fungi    Chlorophyta Dinoflagellata 
#0.2326875      0.4266730      0.5348663      0.6413221      0.7376608 

#$`2_umiujaq_7_ice_shelves`
#Metazoa     Ochrophyta    Chlorophyta Dinoflagellata 
#0.3020749      0.4754613      0.6154493      0.7183333 

#$`3_cambridge_bay_4_bylot_island`
#Chlorophyta     Metazoa 
#0.4322992   0.7189061 

#$`3_cambridge_bay_5_resolute`
#Metazoa  Ochrophyta       Fungi Chlorophyta 
#0.3311679   0.4891919   0.6230541   0.7407043 

#$`3_cambridge_bay_6_ellesmere_island`
#Metazoa  Ochrophyta       Fungi Chlorophyta 
#0.3522202   0.5014681   0.6328242   0.7503479 

#$`3_cambridge_bay_7_ice_shelves`
#Metazoa     Ochrophyta    Chlorophyta          Fungi Dinoflagellata 
#0.2430345      0.4003546      0.5433527      0.6862829      0.7675477 

#$`4_bylot_island_5_resolute`
#Chlorophyta     Metazoa 
#0.4469602   0.8273426 

#$`4_bylot_island_6_ellesmere_island`
#Chlorophyta     Metazoa  Ochrophyta 
#0.4680949   0.6444478   0.8044647 

#$`4_bylot_island_7_ice_shelves`
#Chlorophyta     Metazoa 
#0.4455539   0.7346752 

#$`5_resolute_6_ellesmere_island`
#Metazoa Ochrophyta      Fungi 
#0.4686620  0.6213855  0.7402441 

#$`5_resolute_7_ice_shelves`
#Metazoa Ochrophyta      Fungi 
#0.3373009  0.5989800  0.7129690 

#$`6_ellesmere_island_7_ice_shelves`
#Metazoa  Ochrophyta       Fungi Chlorophyta 
#0.3700136   0.5132044   0.6420063   0.7500815 

#water body simper####
simp_18S_water_body <- simper(ps_18S_ave_RA_Division_count, simper_meta_18S$water_body, permutations = 9999)
simp_18S_water_body

#SIMPER 18s WATER BODY OUTPUT####
#cumulative contributions of most influential species:
#cumulative contributions of most influential species:

#$pond_lake
#Metazoa    Chlorophyta     Ochrophyta Dinoflagellata 
#0.3183691      0.4770312      0.6211863      0.7216841 

#$pond_rock_pool
#Ochrophyta     Metazoa Chlorophyta 
#0.2549790   0.4927117   0.7282120 

#$pond_meltwater_pond
#Metazoa    Chlorophyta     Ochrophyta Dinoflagellata 
#0.3129025      0.4960152      0.6550098      0.7670179 

#$lake_rock_pool
#Metazoa Chlorophyta  Ochrophyta 
#0.3961448   0.6442062   0.8499332 

#$lake_meltwater_pond
#Metazoa  Ochrophyta       Fungi Chlorophyta Apicomplexa 
#0.2454337   0.4690235   0.5822616   0.6823597   0.7609270 

#$rock_pool_meltwater_pond
#Metazoa Chlorophyta 
#0.4215120   0.7244267   

#SIMPER PROTISTS ONLY####
#Finally let's look at simper when we look at just the protist community!
PROTIST_func_ave_RA
unique(tax_table(PROTIST_func_ave_RA)[,"Division"] ) #this tells you that there are 28 unique phyla
PROTIST_func_ave_RA_Division <- tax_glom(PROTIST_func_ave_RA, "Division", NArm = TRUE)
PROTIST_func_ave_RA_Division #this shows us all our ASVs have been sorted into their respective 27 taxa now as there is no metazoa or fungi
unique(tax_table(PROTIST_func_ave_RA_Division)[,"Division"] )

#Now we replace the ASVIDs with the Division names!
taxa_names(PROTIST_func_ave_RA_Division) <- tax_table(PROTIST_func_ave_RA_Division)[,"Division"]
taxa_names(PROTIST_func_ave_RA_Division)

#Make a count table
PROTIST_func_ave_RA_Division_count <- abundances(PROTIST_func_ave_RA_Division) #this utilises the compositional function in microbiome package to create relative abundances of our data
PROTIST_func_ave_RA_Division_count <- as.data.frame(PROTIST_func_ave_RA_Division_count)
#This needs to be transposed so the samples are rows for vegan
PROTIST_func_ave_RA_Division_count <- t(PROTIST_func_ave_RA_Division_count)
rownames(PROTIST_func_ave_RA_Division_count) #good are rows are now our samples

#region simper####
simp_PROTIST_region <- simper(PROTIST_func_ave_RA_Division_count, simper_meta_18S$region, permutations = 9999)
simp_PROTIST_region

#Simper Region Protist Output####
#cumulative contributions of most influential species:
  
#$`1_kuujjuarapik_2_umiujaq`
#Dinoflagellata     Ochrophyta    Chlorophyta       Cercozoa 
#0.2329278      0.4373675      0.6362278      0.7163703 

#$`1_kuujjuarapik_3_cambridge_bay`
#Ochrophyta    Chlorophyta Dinoflagellata 
#0.2685344      0.5191820      0.7237499 

#$`1_kuujjuarapik_4_bylot_island`
#Chlorophyta Dinoflagellata     Ochrophyta 
#0.5466797      0.6888216      0.8133083 

#$`1_kuujjuarapik_5_resolute`
#Chlorophyta     Ochrophyta Dinoflagellata 
#0.2779098      0.4991699      0.7163202 

#$`1_kuujjuarapik_6_ellesmere_island`
#Ochrophyta Dinoflagellata    Chlorophyta 
#0.3120377      0.5150764      0.7105510 

#$`1_kuujjuarapik_7_ice_shelves`
#Ochrophyta    Chlorophyta Dinoflagellata 
#0.2996928      0.5541097      0.7639845 

#$`2_umiujaq_3_cambridge_bay`
#Ochrophyta    Chlorophyta Dinoflagellata       Cercozoa     Ciliophora 
#0.2478989      0.4413901      0.5732890      0.6791884      0.7775247 

#$`2_umiujaq_4_bylot_island`
#Chlorophyta     Ochrophyta Dinoflagellata 
#0.5785556      0.6850536      0.7723226 

#$`2_umiujaq_5_resolute`
#Chlorophyta     Ochrophyta Dinoflagellata       Cercozoa     Ciliophora 
#0.2219729      0.4390475      0.5833820      0.6967765      0.7966734 

#$`2_umiujaq_6_ellesmere_island`
#Ochrophyta    Chlorophyta Dinoflagellata       Cercozoa     Ciliophora 
#0.3030077      0.4511913      0.5873111      0.6878022      0.7808750 

#$`2_umiujaq_7_ice_shelves`
#Ochrophyta    Chlorophyta Dinoflagellata       Cercozoa 
#0.2824676      0.4853846      0.6351809      0.7264614 

#$`3_cambridge_bay_4_bylot_island`
#Chlorophyta  Ochrophyta 
#0.6514222   0.8307376 

#$`3_cambridge_bay_5_resolute`
#Ochrophyta Chlorophyta Apicomplexa  Ciliophora 
#0.3073100   0.4943524   0.6492844   0.7499650 

#$`3_cambridge_bay_6_ellesmere_island`
#Ochrophyta Chlorophyta    Cercozoa Apicomplexa 
#0.2843637   0.4980892   0.6496022   0.7567240 

#$`3_cambridge_bay_7_ice_shelves`
#Ochrophyta    Chlorophyta Dinoflagellata    Apicomplexa 
#0.2575652      0.4617844      0.5993062      0.7321606 

#$`4_bylot_island_5_resolute`
#Chlorophyta 
#0.7584011 

#$`4_bylot_island_6_ellesmere_island`
#Chlorophyta  Ochrophyta 
#0.6371755   0.8515050 

#$`4_bylot_island_7_ice_shelves`
#Chlorophyta  Ochrophyta 
#0.6546819   0.8566092 

#$`5_resolute_6_ellesmere_island`
#Ochrophyta Chlorophyta    Cercozoa Apicomplexa 
#0.3563272   0.5438487   0.6889064   0.7978567 

#$`5_resolute_7_ice_shelves`
#Ochrophyta Dinoflagellata     Ciliophora 
#0.4769485      0.6326525      0.7425826 

#$`6_ellesmere_island_7_ice_shelves`
#Ochrophyta Chlorophyta    Cercozoa Apicomplexa 
#0.2842677   0.4952107   0.6112641   0.7110600 


#water body simper####
simp_PROTIST_water_body <- simper(PROTIST_func_ave_RA_Division_count, simper_meta_18S$water_body, permutations = 9999)
simp_PROTIST_water_body

#Simper Water Body PROTIST Output####
#cumulative contributions of most influential species:
  
#$pond_lake
#Ochrophyta    Chlorophyta Dinoflagellata    Apicomplexa 
#0.2503344      0.4921473      0.6369273      0.7234150 

#$pond_rock_pool
#Ochrophyta Chlorophyta 
#0.3592473   0.7020500 

#$pond_meltwater_pond
#Chlorophyta     Ochrophyta Dinoflagellata 
#0.2736529      0.5434116      0.7082976 

#$lake_rock_pool
#Chlorophyta  Ochrophyta 
#0.4351996   0.8030061 

#$lake_meltwater_pond
#Ochrophyta Chlorophyta Apicomplexa  Ciliophora 
#0.3537502   0.5003099   0.6249109   0.7448408 

#$rock_pool_meltwater_pond
#Chlorophyta  Ochrophyta 
#0.5448060   0.8217888 


#To display these results visually, I often use a boxplot to show the differences in distribution of these identified, statistically significant species between my groupings -> https://jkzorz.github.io/2019/07/02/boxplots.html

#END OF SCRIPT####

