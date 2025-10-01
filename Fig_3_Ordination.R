#V4 16S and V9 18S ordination 
#Author: Patrick M. Hooper
#Created: 30/08/22

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

#OPTIONAL: INSTALL AND LOAD PALETTES##
library("chroma")
library("viridis")

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
meta_18S <- read.table("~/18S_metadata_averaged.txt") #NOTE THIS IS THE METADATA WITH THE READ COUNTS WITHOUT METAZOA AND PLANTS
meta_18S
rownames(meta_18S)
nrow(meta_18S) #32

#6. LOAD YOUR ASV COUNT TABLE .TXT FILE
#I made this in excel from the functional mothership to make sure that ASVIDs were all matching
#No KJB17bii in this file
#This needs to not have a column header on the ASVID column
count_18S <- read.table("~/18S_asv_count_average.txt", header = TRUE) 

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

#Load the functional taxonomy table for the 18S data
func_tax_18S <- read.table("~/18S_functional_taxonomy_table.txt", header = TRUE)
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

#quick check all names are correct
taxa_names(TAX_18S_ave)
taxa_names(OTU_18S_ave)
sample_names(META_18S_ave)
sample_names(OTU_18S_ave)

#4. Loading 16S phyloseq objects
#We are using count tables averaged by water bodies. Biological replicates from water bodies were summed and averaged by the number of replicates. The metadata was also averaged. 

#3. Add averaged metadata
#No KJB17bii in this file
meta_16S <- read.table("~/16S_metadata_average.txt")
head(meta_16S)
rownames(meta_16S)
nrow(meta_16S) #32

#6. LOAD YOUR ASV COUNT TABLE .TXT FILE
#I made this in excel from the functional mothership to make sure that ASVIDs were all matching
#No KJB17bii in this file
#This needs to not have a column header on the ASVID column
count_16S <- read.table("~/FINAL_ps_16S_count_glom_0_01_AVERAGE.txt", header = TRUE) 

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
tax_16S <- read.table("~/FINAL_ps_16S_tax_no_KJ17bii.txt", header = TRUE)
head(tax_16S)

#use this function to add ASVIDS as rownames on the tax table
rownames(tax_16S) <- tax_16S[,1]
rownames(tax_16S)

#Need to remove the column ASVID
tax_16S <- tax_16S %>% select(-ASVID)

#I also had to change it to a matrix
tax_16S <- as.matrix(tax_16S)

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

#remove 0 ASVs
#Check that there are no ASVs that are not present in any of the samples
ps_16S_ave <- prune_taxa(taxa_sums(ps_16S_ave) > 0, ps_16S_ave)
ps_16S_ave #this should leave us with 15250 by 32 samples

#quick check all names are correct
taxa_names(TAX_16S_ave)
taxa_names(OTU_16S_ave)
sample_names(META_16S_ave)
sample_names(OTU_16S_ave)

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

#Figure 3####
set.seed(1442)

                                        
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
                                        
#Bray-Curtis plots                                      
ps_func_18S_ave_RA_BCNMDS <- ordinate(ps_func_18S_ave_RA, distance = "bray", method = "NMDS")
ps_func_18S_ave_RA_BCNMDS

ps_16S_ave_RA_BCNMDS <- ordinate(ps_16S_ave_RA, distance = "bray", method = "NMDS")
ps_16S_ave_RA_BCNMDS 

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

#Anosim calculations####
#check for significant difference between regions with anosim

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


#END OF SCRIPT####

