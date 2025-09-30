#Alpha diversity analysis of V4 16S reads
#Author: Patrick M. Hooper
#Date Created: 16/08/22
#Date Modified: 06/09/22

#1. LOAD PACKAGES####
library(tidyverse); packageVersion("tidyverse")
library(plyr); packageVersion("plyr")
library(dplyr); packageVersion("dplyr")
library(phyloseq); packageVersion("phyloseq")
library(RColorBrewer); packageVersion("RColorBrewer")
library(microbiome); packageVersion("microbiome")
library(vegan); packageVersion("vegan")
library(DESeq2); packageVersion("DESeq2")

##OPTIONAL: LOAD PALETTES
library("chroma")
library("viridis")

#2. Colour palettes and ggplot
#OPTIONAL: set ggplot2 theme to minimal
theme_set(theme_minimal())

#Set standardized plot font sizes
font_size <- theme(axis.text = element_text(size = 18)) + theme(axis.title = element_text(size = 18)) + theme(plot.title = element_text(size = 20)) + theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size = 20)) + theme(plot.subtitle = element_text(size = 18))

#My palette is an adaptation of Pal Tol's color-blind friendly color blind scheme

#SOURCE: Paul Tol: https://personal.sron.nl/~pault/
#Tol_bright <- c('#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')

#More info on palettes here: https://thenode.biologists.com/data-visualization-with-flying-colors/research/

#Set custom colour palette for qualitative data (i moved the colours around to reflect the ecozones...ish)
region_palette <- c("1_kuujjuarapik" = "#228833", "2_umiujaq" = "#CCBB44", "3_cambridge_bay" = "#EE6677","4_bylot_island" = "#AA3377", "5_resolute" = "#66CCEE", "6_ellesmere_island" = "#4477AA", "7_ice_shelves" = "#BBBBBB")
show_col(region_palette)

#For scaled data I will use a two tone color palette, designed to be accessible for most eyesights and available for integration with ggplot2 from my color palette above

#3. Set working directory####
setwd("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS")

#4. Load metadata table####
#This metadata table includes the final version sample names####
meta_16S <- read.table("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/GLOM_V4_16S_METADATA_FORMATTED_UPDATE_17_8_22.txt")

#5. Load phyloseq objects####

#4. LOAD YOUR FINAL 16S PHYLOSEQ TABLE####
ps_16S <- readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FINAL_ps_16S")
ps_16S 

#5. UPDATE METADATA IN PHYLOSEQ OBJECT WITH NEW VARIABLES####
#turn imported meta .txt file into a sample_data dataframe for phyloseq
meta_16S  <- sample_data(meta_16S)
sample_data(ps_16S) <- meta_16S
ps_16S #this should now have 26 metadata variables as we have updated the sample_data in our previous phyloseq object with the new sample names and some updates to the metadata, this metadata has the readcount already added

#Remove sample KJ17bii as we have four replicates for this sample
ps_16S <- prune_samples(sample_names(ps_16S) != "KJ17Bii", ps_16S)
ps_16S #should now have 1520 taxa, 96 samples (4 more than 16S), and 26 sample variables

#Check that there are no ASVs that are not present in any of the samples
ps_16S <- prune_taxa(taxa_sums(ps_16S) > 0, ps_16S)
ps_16S #this should not lose any ASVs

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 15250 taxa and 96 samples ]
#sample_data() Sample Data:       [ 96 samples by 26 sample variables ]
#tax_table()   Taxonomy Table:    [ 15250 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 15250 tips and 15249 internal nodes ]
#refseq()      DNAStringSet:      [ 15250 reference sequences ]


#Summarise sample depth in your phyloseq object across each sample
summary_16S <- as.data.frame(sort(sample_sums(ps_16S)))
summary_16S
summary(summary_16S)

#> summary(summary_16S)
#sort(sample_sums(ps_16S))
#Min.   :  9752           
#1st Qu.:100363           
#Median :125972           
#Mean   :123287           
#3rd Qu.:141073           
#Max.   :263040   

#Total ASV sequences
sum(summary_16S)
#11,835,544

#It would also be useful to know how these samples distribute across the different samples!
#Let's quickly make a histogram of the sample count across each sample in our dataset
hist(sample_sums(ps_16S), main="Histogram: Read Counts", xlab="Total Reads", las=1, breaks=12)


#We are now ready to start our analysis using this phyloseq object

#8. ALPHA DIVERSITY ANALYSIS of 16S ####
#Plotting richness and diversity####
#McMurdie: I know it is tempting to trim noise right away, but many richness estimates are modeled on singletons and doubletons in the abundance data. You need to leave them in the dataset if you want a meaningful estimate.
#Source: https://joey711.github.io/phyloseq/plot_richness-examples.html

#However, DADA2 removes singletons as part of ASV assembly so these can't be included
#Source: https://github.com/benjjneb/dada2/issues/320

#As such, Ben recommends to NOT use richness estimators on ASVs.. https://github.com/benjjneb/dada2/issues/317
#An alternative would be to use Amy Willis' breakaway package: https://adw96.github.io/breakaway/articles/breakaway.html 

#For now I am going to calculate in phyloseq and come back to this later using the ACE index. This index a nonparametric method for estimating the number of species using sample coverage,which is defined as the sum of the probabilities of the observed species. The ACE method divides observed frequencies into abundant and rare groups. The abundant species are those with more than 10 individuals in the sample, and the rare species are those with fewer than 10 individuals. Only the presence or absence information of abundant species is considered in the ACE method because they would be discovered anyway. Therefore, the exact frequencies for the abundant species are not required in the ACE method. On the other hand, the exact frequencies for the rare species are required because the estimation of the number of missing species is based entirely on these rare species. I will use this over Chao which only uses singletons and doubletons and therefore would be less accurate for dada2 data. 

#The uses of Shannon-Weaver and Simpson diversity indices have been recommended to robustly measure microbial diversity, we shall use both:
#Shannon-Weaver = Estimator of species richness and species evenness: more weight on species richness
#Simpson = Estimator of species richness and species evenness: more weight on species evenness
#Reference: Kim et al. 2017, 

#Amy Willis does not recommend rarefying prior to alpha diversity so we are going to use our original phyloseq object for now
#Source: https://www.frontiersin.org/articles/10.3389/fmicb.2019.02407/full 

#Alpha diversity of un-rarefied 16S ASVs ####

#An easy way to estimate richness across multiple diversity indices in phyloseq
#Add latitude and water temperature data from your metadata, this will allow you to do a regression analysis later on
richness_16S <- estimate_richness(ps_16S, measures = c("ACE", "Shannon", "InvSimpson"))%>%
  cbind(sample_data(ps_16S)[,"latitude"])%>%
  cbind(sample_data(ps_16S)[,"water_temp"])
print(richness_16S)

summary(richness_16S)

#PLOT 1: Diversity measures against latitude ####
#1.1 ACE
#Ace plots give a standard error range as part of the calculation
latitude_ACE_16S <- plot_richness(ps_16S, x="latitude", color = "ave_ann_temp", measures=c("ACE"), title = "ACE diversity index of 16S samples plotted by latitude") +  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank()) + geom_point(size=8, alpha=0.8, aes(shape = factor(water_body))) + geom_smooth(method = "lm", alpha = 0.2, colour = "black")  # as the relationship looks linear we can use LM, the grey line indicates the 95% confidence interval
latitude_ACE_16S $layers # useful tool for seeing the layers on your ggplot
latitude_ACE_16S $layers <- latitude_ACE_16S $layers[-1] #removes the small dots that are produced by the phyloseq command so we can use geom_point instead

#replace colour scale with custom color scale to match figures
latitude_ACE_16S <- latitude_ACE_16S + scale_colour_continuous(high = "#EE6677", low = "#4477AA")

#Add linear model
#It's looking like there is a relationship between the diversity indices and latitude. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_16S$ACE~richness_16S$latitude)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.2056, p value = 2.086e-06

latitude_ACE_16S <- latitude_ACE_16S + labs(subtitle = "Adjusted R-squared value: 0.2056, p value <0.001")
latitude_ACE_16S  + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) #check you're happy

#use ggsave to save to your wd as a pdf and jpeg
ggsave("latitude_ACE_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("latitude_ACE_16S.pdf", width = 297, height = 210, units = c("mm"))

#1.2 SHANNON
#Both indexes are used to measure similar concepts of alpha diversity (Simpson's index is less sensitive to the difference in taxa richness than Shannon's index); however, the interpretation is inverse. The lower value of Shannon's index, the lower diversity. The lower value of Simpson's index (range: 0-1), the higher diversity. Since this is quite a non-intuitive scale, the inverse Simpson index is more frequently reported.
latitude_Shannon_16S <- plot_richness(ps_16S, x="latitude", color = "ave_ann_temp", measures=c("Shannon"), title = "Shannon diversity index of 16S samples plotted by latitude") +  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank()) + geom_point(size=8, alpha=0.8, aes(shape = factor(water_body))) + geom_smooth(method = "lm", alpha = 0.2, colour = "black")
latitude_Shannon_16S $layers
latitude_Shannon_16S $layers <- latitude_Shannon_16S $layers[-1]
latitude_Shannon_16S <- latitude_Shannon_16S + scale_colour_continuous(high = "#EE6677", low = "#4477AA")

#Add linear model
#It's looking like there is a relationship between the diversity indices and latitude. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_16S$Shannon~richness_16S$latitude)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.2152, p value = 1.153e-06

latitude_Shannon_16S <- latitude_Shannon_16S + labs(subtitle = "Adjusted R-squared value: 0.2152 , p value <0.001")
latitude_Shannon_16S  + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("latitude_Shannon_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("latitude_Shannon_16S.pdf", width = 297, height = 210, units = c("mm"))

#1.2 INVERSE SIMPSON
#Both indexes are used to measure similar concepts of alpha diversity (Simpson's index is less sensitive to the difference in taxa richness than Shannon's index); however, the interpretation is inverse. The lower value of Shannon's index, the lower diversity. The lower value of Simpson's index (range: 0-1), the higher diversity. Since this is quite a non-intuitive scale, the inverse Simpson index is more frequently reported.
latitude_InvSimpson_16S <- plot_richness(ps_16S, x="latitude", color = "ave_ann_temp", measures=c("InvSimpson"), title = "InvSimpson diversity index of 16S samples plotted by latitude") +  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank()) + geom_point(size=8, alpha=0.8, aes(shape = factor(water_body))) + geom_smooth(method = "lm", alpha = 0.2, colour = "black")
latitude_InvSimpson_16S $layers
latitude_InvSimpson_16S $layers <- latitude_InvSimpson_16S $layers[-1]
latitude_InvSimpson_16S <- latitude_InvSimpson_16S + scale_colour_continuous(high = "#EE6677", low = "#4477AA")

#Add linear model
#It's looking like there is a relationship between the diversity indices and latitude. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_16S$InvSimpson~richness_16S$latitude)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.1626, p value = 2.75e-05

latitude_InvSimpson_16S <- latitude_InvSimpson_16S + labs(subtitle = "Adjusted R-squared value: 0.1626, p value <0.001")
latitude_InvSimpson_16S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("latitude_InvSimpson_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("latitude_InvSimpson_16S.pdf", width = 297, height = 210, units = c("mm"))

#PLOT 2. WATER TEMPERATURE####
#Note Resolute Shore Pond does not have water temperature data
#2.1 ACE
#Ace plots give a standard error range as part of the calculation

water_temp_ACE_16S<-plot_richness(ps_16S, x="water_temp", measures=c("ACE"), title = "ACE diversity index of 16S samples plotted against water temperature")+ geom_point(size=8, alpha=0.8, aes(colour = factor(region), shape=factor(water_body))) + theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())+
  geom_smooth(colour = "black", method="lm", size=0.5)
water_temp_ACE_16S $layers
water_temp_ACE_16S $layers <- water_temp_ACE_16S $layers[-1]
water_temp_ACE_16S + scale_color_manual(values = region_palette)

#Add linear model
#It's looking like there is a relationship between the diversity indices and water temperature. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_16S$ACE~richness_16S$water_temp)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.2124, p value=1.997e-06

water_temp_ACE_16S <- water_temp_ACE_16S + labs(subtitle = "Adjusted R-squared value: 0.2124, p value <0.001")+ scale_color_manual(values = region_palette)
water_temp_ACE_16S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("water_temp_ACE_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("water_temp_ACE_16S.pdf", width = 297, height = 210, units = c("mm"))

water_temp_Shannon_16S<-plot_richness(ps_16S, x="water_temp", measures=c("Shannon"), title = "Shannon diversity index of 16S samples plotted against water temperature")+ geom_point(size=8, alpha=0.8, aes(colour = factor(region), shape=factor(water_body))) + theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())+
  geom_smooth(colour = "black", method="lm", size=0.5)
water_temp_Shannon_16S $layers
water_temp_Shannon_16S $layers <- water_temp_Shannon_16S $layers[-1]
water_temp_Shannon_16S

#Add linear model
#It's looking like there is a relationship between the diversity indices and water temperature. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_16S$Shannon~richness_16S$water_temp)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.1711, p value=2.244e-05

water_temp_Shannon_16S <- water_temp_Shannon_16S + labs(subtitle = "Adjusted R-squared value: 0.1711, p value <0.001")+ scale_color_manual(values = region_palette)
water_temp_Shannon_16S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("water_temp_Shannon_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("water_temp_Shannon_16S.pdf", width = 297, height = 210, units = c("mm"))

water_temp_InvSimpson_16S<-plot_richness(ps_16S, x="water_temp", measures=c("InvSimpson"), title = "Inverse Simpson diversity index of 16S samples plotted against water temperature")+ geom_point(size=8, alpha=0.8, aes(colour = factor(region), shape=factor(water_body))) + theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())+
  geom_smooth(colour = "black", method="lm", size=0.5)
water_temp_InvSimpson_16S $layers
water_temp_InvSimpson_16S $layers <- water_temp_InvSimpson_16S $layers[-1]
water_temp_InvSimpson_16S+ scale_color_manual(values = region_palette)

#Add linear model
#It's looking like there is a relationship between the diversity indices and water temperature. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_16S$InvSimpson~richness_16S$water_temp)
summary(lm)
s#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.07792, p value=0.003896

water_temp_InvSimpson_16S <- water_temp_InvSimpson_16S + labs(subtitle = "Adjusted R-squared value: 0.07792, p value = 0.00390")+ scale_color_manual(values = region_palette)
water_temp_InvSimpson_16S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("water_temp_InvSimpson_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("water_temp_InvSimpson_16S.pdf", width = 297, height = 210, units = c("mm"))

#PLOT 3. Comparing alpha diversity variation within pond samples
#This is the names of the samples ordered by latitude
summed_order = c("KJ1","KJ2","KJ3","KJ4","KJ5","KJ6","KJ7","KJ8","UM1","UM2","UM3","UM4","UM5","UM6","UM7","UM8","CB1","CB2","CB3","CB4","CB5","CB6","CB7","CB8","BY1","BY2","RE1","RE2","WH","AP","MKIS","WHIS")

ponds_ACE_16S <- plot_richness(ps_16S, x="new_pond", color ="region", shape = "water_body", measures=c("ACE"), title = "ACE diversity index of 16S samples grouped by sample sites")+
  geom_point(size=8, alpha=0.8)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())  +
  geom_boxplot(alpha=0.6)
ponds_ACE_16S $layers
ponds_ACE_16S $layers <- ponds_ACE_16S $layers[-1]
ponds_ACE_16S <- ponds_ACE_16S + scale_x_discrete(limits=summed_order) + scale_color_manual(values = region_palette)
ponds_ACE_16S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("ponds_ACE_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("ponds_ACE_16S.pdf", width = 297, height = 210, units = c("mm"))

ponds_Shannon_16S <- plot_richness(ps_16S, x="new_pond", color ="region", shape = "water_body", measures=c("Shannon"),  title = "Shannon diversity index of 16S samples grouped by sample sites")+
  geom_point(size=8, alpha=0.8)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())  +
  geom_boxplot(alpha=0.6)
ponds_Shannon_16S $layers
ponds_Shannon_16S $layers <- ponds_Shannon_16S $layers[-1]
ponds_Shannon_16S <- ponds_Shannon_16S + scale_x_discrete(limits=summed_order) + scale_color_manual(values = region_palette)
ponds_Shannon_16S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("ponds_Shannon_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("ponds_Shannon_16S.pdf", width = 297, height = 210, units = c("mm"))

ponds_InvSimpson_16S <- plot_richness(ps_16S, x="new_pond", color ="region", shape = "water_body", measures=c("InvSimpson"),  title = "Inverse Simpson diversity index of 16S samples diversity grouped by sample sites")+
  geom_point(size=8, alpha=0.8)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())  +
  geom_boxplot(alpha=0.6)
ponds_InvSimpson_16S $layers
ponds_InvSimpson_16S $layers <- ponds_InvSimpson_16S $layers[-1]
ponds_InvSimpson_16S <- ponds_InvSimpson_16S + scale_x_discrete(limits=summed_order) + scale_color_manual(values = region_palette)
ponds_InvSimpson_16S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("ponds_InvSimpson_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("ponds_InvSimpson_16S.pdf", width = 297, height = 210, units = c("mm"))

#PLOT 4. Alpha diversity by regions
region_ACE_16S <- plot_richness(ps_16S, x="region", color ="region", measures=c("ACE"), title = "ACE diversity index of 16S samples grouped by region")+
  geom_jitter(size=8, alpha=0.8, aes(shape  = water_body), width = 0.2)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())
region_ACE_16S $layers
region_ACE_16S $layers <- region_ACE_16S $layers[-1]
region_ACE_16S <- region_ACE_16S + scale_color_manual(values = region_palette) + geom_boxplot(width=0.1, alpha=0.8)
region_ACE_16S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("region_ACE_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("region_ACE_16S.pdf", width = 297, height = 210, units = c("mm"))

region_Shannon_16S <- plot_richness(ps_16S, x="region", color ="region", measures=c("Shannon"), title = "Shannon diversity index of 16S samples grouped by region")+
  geom_jitter(size=8, alpha=0.8, aes(shape  = water_body), width = 0.2)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())
region_Shannon_16S $layers
region_Shannon_16S $layers <- region_Shannon_16S $layers[-1]
region_Shannon_16S <- region_Shannon_16S + scale_color_manual(values = region_palette) + geom_boxplot(width=0.1, alpha=0.8)
region_Shannon_16S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("region_Shannon_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("region_Shannon_16S.pdf", width = 297, height = 210, units = c("mm"))

region_InvSimpson_16S <- plot_richness(ps_16S, x="region", color ="region", measures=c("InvSimpson"), title = "Inverse Simpson diversity index of 16S samples grouped by region")+
  geom_jitter(size=8, alpha=0.8, aes(shape  = water_body), width = 0.2)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())
region_InvSimpson_16S $layers
region_InvSimpson_16S $layers <- region_InvSimpson_16S $layers[-1]
region_InvSimpson_16S <- region_InvSimpson_16S + scale_color_manual(values = region_palette) + geom_boxplot(width=0.1, alpha=0.8)
region_InvSimpson_16S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("region_InvSimpson_16S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("region_InvSimpson_16S.pdf", width = 297, height = 210, units = c("mm"))

#calculating richness and diversity####
#Null hypothesis = no difference in richness or evenness between sample regions.

#We will use the estimate_richness table we made earlier with the latitude and water temperature
richness_16S 

#We will add the region and new_ID metadata
richness_16S <- richness_16S %>% cbind(sample_data(ps_16S)[,"region"])%>% cbind(sample_data(ps_16S)[,"new_ID"])
richness_16S

#Save this output as a .txt file
write.table(richness_16S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/final_results/richness_estimates_16S.txt")

#Check skew of data (normality)
hist(richness_16S$ACE)
hist(richness_16S$Shannon)
hist(richness_16S$InvSimpson)

#Test for normality (If p >0.05, data is normal)
shapiro.test(richness_16S$ACE) #p = 0.2174
shapiro.test(richness_16S$Shannon) #p = 3.541e-06
shapiro.test(richness_16S$InvSimpson) #p = 8.449e-09
#This shows us that the ACE data is normal but the evenness values for Shannon and Simpson are not normal

#As such we will use the non-parametric test for all diversity measures for continuity 

kruskal.test(richness_16S$ACE ~ richness_16S$region, data = richness_16S)
#data:  richness_16S$ACE by richness_16S$region
#Kruskal-Wallis chi-squared = 42.589, df = 6,
#p-value = 1.406e-07

ACE_16S_wilcoxon_summary <- pairwise.wilcox.test(richness_16S$ACE, richness_16S$region, p.adjust.method = "none")
ACE_16S_wilcoxon_summary

Pairwise comparisons using Wilcoxon rank sum test with continuity correction 

#data:  richness_16S$ACE and richness_16S$region 

#                   1_kuujjuarapik 2_umiujaq
#2_umiujaq          0.00062        -        
#3_cambridge_bay    3.0e-06        0.91055  
#4_bylot_island     0.00082        0.59527  
#5_resolute         0.00758        0.63095  
#6_ellesmere_island 0.01831        0.29623  
#7_ice_shelves      0.00021        3.4e-06  

#                   3_cambridge_bay 4_bylot_island
#2_umiujaq          -               -             
#3_cambridge_bay    -               -             
# 4_bylot_island    0.27619         -             
#5_resolute         0.59527         0.13203       
#6_ellesmere_island 0.02464         0.00433       
#7_ice_shelves      3.4e-06         0.00216      

#                   5_resolute 6_ellesmere_island
#2_umiujaq          -          -                 
#3_cambridge_bay    -          -                 
#4_bylot_island     -          -                 
#5_resolute         -          -                 
#6_ellesmere_island 0.48485    -                 
#7_ice_shelves      0.00216    0.00216  

kruskal.test(richness_16S$Shannon ~ richness_16S$region, data = richness_16S)
#data:  richness_16S$Shannon by richness_16S$region
#Kruskal-Wallis chi-squared = 40.633, df = 6, p-value = 3.42e-07

Shannon_16S_wilcoxon_summary <- pairwise.wilcox.test(richness_16S$Shannon, richness_16S$region, p.adjust.method = "none")
Shannon_16S_wilcoxon_summary

#data:  richness_16S$Shannon and richness_16S$region 

#                   1_kuujjuarapik 2_umiujaq
#2_umiujaq          0.00117        -        
#3_cambridge_bay    5.4e-06        0.87824  
#4_bylot_island     0.00010        0.02864  
#5_resolute         0.01529        0.93966  
#6_ellesmere_island 0.00906        0.85975  
#7_ice_shelves      3.4e-06        0.00015  
#                   3_cambridge_bay 4_bylot_island
#2_umiujaq          -               -             
#3_cambridge_bay    -               -               #4_bylot_island     0.00042         -             
#5_resolute         0.93966         0.13203       
#6_ellesmere_island 0.66743         0.01515       
#7_ice_shelves      1.3e-05         0.24026       
#                   5_resolute 6_ellesmere_island
#2_umiujaq          -          -                 
#3_cambridge_bay    -          -                 
#4_bylot_island     -          -                 
#5_resolute         -          -                 
#6_ellesmere_island 1.00000    -                 
#7_ice_shelves      0.00216    0.00216 
  

kruskal.test(richness_16S$InvSimpson ~ richness_16S$region, data = richness_16S)
#data:  richness_16S$InvSimpson by richness_16S$region
#Kruskal-Wallis chi-squared = 29.964, df = 6, p-value = 3.994e-05

InvSimpson_16S_wilcoxon_summary <- pairwise.wilcox.test(richness_16S$InvSimpson, richness_16S$region, p.adjust.method = "none")
InvSimpson_16S_wilcoxon_summary
#data:  richness_16S$InvSimpson and richness_16S$region 

#                   1_kuujjuarapik 2_umiujaq
#2_umiujaq          0.00823        -        
#3_cambridge_bay    0.00085        1.00000  
#4_bylot_island     0.00022        0.02464  
#5_resolute         0.05699        0.89956  
#6_ellesmere_island 0.02864        0.89956  
#7_ice_shelves      0.00030        0.00906  
#                   3_cambridge_bay 4_bylot_island
#2_umiujaq          -               -             
#3_cambridge_bay    -               -             
#4_bylot_island     0.00168         -             
#5_resolute         0.82024         0.09307       
#6_ellesmere_island 0.70470         0.04113       
#7_ice_shelves      0.00057         0.81818       
#                   5_resolute 6_ellesmere_island
#2_umiujaq          -          -                 
#3_cambridge_bay    -          -                 
#4_bylot_island     -          -                 
#5_resolute         -          -                 
#6_ellesmere_island 0.81818    -                 
# _ice_shelves      0.00216    0.01515    


