#Alpha diversity analysis of V9 18S reads
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
#This metadata contains information on all samples except the surplus KJ17bii
meta_18S <- read.table("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/GLOM_V9_18S_METADATA_FORMATTED_UPDATE_17_8_22.txt")

#5. Load phyloseq objects####
#4. LOAD YOUR FINAL 18S PHYLOSEQ TABLE WITHOUT METAZOA OR PLANTS####
ps_18S <- readRDS(file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/saveRDS/FINAL_ps_18S_no_met_plants")
ps_18S

#Remove sample KJ17bii as we have four replicates for this sample
ps_18S <- prune_samples(sample_names(ps_18S) != "KJ17Bii", ps_18S)
ps_18S #should now have 92 samples

#add the read count to ps_18S
read_count <- readcount(ps_18S)
sample_data(ps_18S)$read_count <- read_count
head(sample_data(ps_18S))

#5. UPDATE METADATA IN PHYLOSEQ OBJECT WITH NEW VARIABLES####
#turn imported meta .txt file into a sample_data dataframe for phyloseq
meta_18S  <- sample_data(meta_18S)

#Adds the new sample IDs to the phyloseq object
sample_data(ps_18S) <- meta_18S

ps_18S #this should now have 8179 taxa, 92 samples, and 25 sample variables


#Summarise sample depth in your phyloseq object across each sample
summary_18S <- as.data.frame(sort(sample_sums(ps_18S)))
summary_18S
summary(summary_18S)

#> summary(summary_18S)
#sort(sample_sums(ps_18S))
#Min.   : 5842            
#1st Qu.:19092            
#Median :26401            
#Mean   :29450            
#3rd Qu.:38102            
#Max.   :68364  

#Total ASV sequences
sum(summary_18S)
#2,709,404

#It would also be useful to know how these samples distribute across the different samples!
#Let's quickly make a histogram of the sample count across each sample in our dataset
hist(sample_sums(ps_18S), main="Histogram: Read Counts", xlab="Total Reads", las=1, breaks=12)


#We are now ready to start our analysis using this phyloseq object

#8. ALPHA DIVERSITY ANALYSIS of 18S ####
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

#Alpha diversity of un-rarefied 18S ASVs ####

#An easy way to estimate richness across multiple diversity indices in phyloseq
#Add latitude and water temperature data from your metadata, this will allow you to do a regression analysis later on
richness_18S <- estimate_richness(ps_18S, measures = c("ACE", "Shannon", "InvSimpson"))%>%
  cbind(sample_data(ps_18S)[,"latitude"])%>%
  cbind(sample_data(ps_18S)[,"water_temp"])
print(richness_18S)


#PLOT 1: Diversity measures against latitude ####
#1.1 ACE
#Ace plots give a standard error range as part of the calculation
latitude_ACE_18S <- plot_richness(ps_18S, x="latitude", color = "ave_ann_temp", measures=c("ACE"), title = "ACE diversity index of 18S samples plotted by latitude") +  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank()) + geom_point(size=8, alpha=0.8, aes(shape = factor(water_body))) + geom_smooth(method = "lm", alpha = 0.2, colour = "black")  # as the relationship looks linear we can use LM, the grey line indicates the 95% confidence interval
latitude_ACE_18S $layers # useful tool for seeing the layers on your ggplot
latitude_ACE_18S $layers <- latitude_ACE_18S $layers[-1] #removes the small dots that are produced by the phyloseq command so we can use geom_point instead

#replace colour scale with custom color scale to match figures
latitude_ACE_18S <- latitude_ACE_18S + scale_colour_continuous(high = "#EE6677", low = "#4477AA")

#Add linear model
#It's looking like there is a relationship between the diversity indices and latitude. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_18S$ACE~richness_18S$latitude)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.3855, p value«<0.001

latitude_ACE_18S <- latitude_ACE_18S + labs(subtitle = "Adjusted R-squared value: 0.3855, p value <0.001")
latitude_ACE_18S#check you're happy
latitude_ACE_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("latitude_ACE_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("latitude_ACE_18S.pdf", width = 297, height = 210, units = c("mm"))

#1.2 SHANNON
#Both indexes are used to measure similar concepts of alpha diversity (Simpson's index is less sensitive to the difference in taxa richness than Shannon's index); however, the interpretation is inverse. The lower value of Shannon's index, the lower diversity. The lower value of Simpson's index (range: 0-1), the higher diversity. Since this is quite a non-intuitive scale, the inverse Simpson index is more frequently reported.
latitude_Shannon_18S <- plot_richness(ps_18S, x="latitude", color = "ave_ann_temp", measures=c("Shannon"), title = "Shannon diversity index of 18S samples plotted by latitude") +  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank()) + geom_point(size=8, alpha=0.8, aes(shape = factor(water_body))) + geom_smooth(method = "lm", alpha = 0.2, colour = "black")
latitude_Shannon_18S $layers
latitude_Shannon_18S $layers <- latitude_Shannon_18S $layers[-1]
latitude_Shannon_18S <- latitude_Shannon_18S + scale_colour_continuous(high = "#EE6677", low = "#4477AA")

#Add linear model
#It's looking like there is a relationship between the diversity indices and latitude. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_18S$Shannon~richness_18S$latitude)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.1726, p value«<0.001

latitude_Shannon_18S <- latitude_Shannon_18S + labs(subtitle = "Adjusted R-squared value: 0.1726, p value <0.001")
latitude_Shannon_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("latitude_Shannon_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("latitude_Shannon_18S.pdf", width = 297, height = 210, units = c("mm"))

#1.2 INVERSE SIMPSON
#Both indexes are used to measure similar concepts of alpha diversity (Simpson's index is less sensitive to the difference in taxa richness than Shannon's index); however, the interpretation is inverse. The lower value of Shannon's index, the lower diversity. The lower value of Simpson's index (range: 0-1), the higher diversity. Since this is quite a non-intuitive scale, the inverse Simpson index is more frequently reported.
latitude_InvSimpson_18S <- plot_richness(ps_18S, x="latitude", color = "ave_ann_temp", measures=c("InvSimpson"), title = "InvSimpson diversity index of 18S samples plotted by latitude") +  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank()) + geom_point(size=8, alpha=0.8, aes(shape = factor(water_body))) + geom_smooth(method = "lm", alpha = 0.2, colour = "black")
latitude_InvSimpson_18S $layers
latitude_InvSimpson_18S $layers <- latitude_InvSimpson_18S $layers[-1]
latitude_InvSimpson_18S <- latitude_InvSimpson_18S + scale_colour_continuous(high = "#EE6677", low = "#4477AA")

#Add linear model
#It's looking like there is a relationship between the diversity indices and latitude. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_18S$InvSimpson~richness_18S$latitude)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.05275, p value«<0.01567

latitude_InvSimpson_18S <- latitude_InvSimpson_18S + labs(subtitle = "Adjusted R-squared value: 0.05275, p value = 0.0157")
latitude_InvSimpson_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("latitude_InvSimpson_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("latitude_InvSimpson_18S.pdf", width = 297, height = 210, units = c("mm"))

#PLOT 2. WATER TEMPERATURE####
#Note Resolute Shore Pond does not have water temperature data
#2.1 ACE
#Ace plots give a standard error range as part of the calculation

water_temp_ACE_18S<-plot_richness(ps_18S, x="water_temp", measures=c("ACE"), title = "ACE diversity index of 18S samples plotted against water temperature")+ geom_point(size=8, alpha=0.8, aes(colour = factor(region), shape=factor(water_body))) + theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())+
  geom_smooth(colour = "black", method="lm", size=0.5)
water_temp_ACE_18S $layers
water_temp_ACE_18S $layers <- water_temp_ACE_18S $layers[-1]
water_temp_ACE_18S + scale_color_manual(values = region_palette)

#Add linear model
#It's looking like there is a relationship between the diversity indices and water temperature. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_18S$ACE~richness_18S$water_temp)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.4156, p value=5.579e-12

water_temp_ACE_18S <- water_temp_ACE_18S + labs(subtitle = "Adjusted R-squared value: 0.4156, p value <0.001")+ scale_color_manual(values = region_palette)
water_temp_ACE_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("water_temp_ACE_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("water_temp_ACE_18S.pdf", width = 297, height = 210, units = c("mm"))

water_temp_Shannon_18S<-plot_richness(ps_18S, x="water_temp", measures=c("Shannon"), title = "Shannon diversity index of 18S samples plotted against water temperature")+ geom_point(size=8, alpha=0.8, aes(colour = factor(region), shape=factor(water_body))) + theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())+
  geom_smooth(colour = "black", method="lm", size=0.5)
water_temp_Shannon_18S $layers
water_temp_Shannon_18S $layers <- water_temp_Shannon_18S $layers[-1]
water_temp_Shannon_18S

#Add linear model
#It's looking like there is a relationship between the diversity indices and water temperature. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_18S$Shannon~richness_18S$water_temp)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.1922, p value=1.032e-05

water_temp_Shannon_18S <- water_temp_Shannon_18S + labs(subtitle = "Adjusted R-squared value: 0.1922, p value <0.001 ")+ scale_color_manual(values = region_palette)
water_temp_Shannon_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("water_temp_Shannon_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("water_temp_Shannon_18S.pdf", width = 297, height = 210, units = c("mm"))

water_temp_InvSimpson_18S<-plot_richness(ps_18S, x="water_temp", measures=c("InvSimpson"), title = "Inverse Simpson diversity index of 18S samples plotted against water temperature")+ geom_point(size=8, alpha=0.8, aes(colour = factor(region), shape=factor(water_body))) + theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())+
  geom_smooth(colour = "black", method="lm", size=0.5)
water_temp_InvSimpson_18S $layers
water_temp_InvSimpson_18S $layers <- water_temp_InvSimpson_18S $layers[-1]
water_temp_InvSimpson_18S+ scale_color_manual(values = region_palette)

#Add linear model
#It's looking like there is a relationship between the diversity indices and water temperature. To statistically support this relationship, perform a linear regression in R to get the R2 and p values:
lm = lm(richness_18S$InvSimpson~richness_18S$water_temp)
summary(lm)
#The output of this shows the relationship is highly statistically significant (Adjusted R-squared value: 0.0973, p value=0.001704

water_temp_InvSimpson_18S <- water_temp_InvSimpson_18S + labs(subtitle = "Adjusted R-squared value: 0.0973, p value = 0.00170")+ scale_color_manual(values = region_palette)
water_temp_InvSimpson_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("water_temp_InvSimpson_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("water_temp_InvSimpson_18S.pdf", width = 297, height = 210, units = c("mm"))

#PLOT 3. Comparing alpha diversity variation within pond samples
#This is the names of the samples ordered by latitude
summed_order = c("KJ1","KJ2","KJ3","KJ4","KJ5","KJ6","KJ7","KJ8","UM1","UM2","UM3","UM4","UM5","UM6","UM7","UM8","CB1","CB2","CB3","CB4","CB5","CB6","CB7","CB8","BY1","BY2","RE1","RE2","WH","AP","MKIS","WHIS")

ponds_ACE_18S <- plot_richness(ps_18S, x="new_pond", color ="region", shape = "water_body", measures=c("ACE"), title = "ACE diversity index of 18S samples grouped by sample sites")+
  geom_point(size=8, alpha=0.8)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())  +
  geom_boxplot(alpha=0.6)
ponds_ACE_18S $layers
ponds_ACE_18S $layers <- ponds_ACE_18S $layers[-1]
ponds_ACE_18S <- ponds_ACE_18S + scale_x_discrete(limits=summed_order) + scale_color_manual(values = region_palette)
ponds_ACE_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("ponds_ACE_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("ponds_ACE_18S.pdf", width = 297, height = 210, units = c("mm"))

ponds_Shannon_18S <- plot_richness(ps_18S, x="new_pond", color ="region", shape = "water_body", measures=c("Shannon"),  title = "Shannon diversity index of 18S samples grouped by sample sites")+
  geom_point(size=8, alpha=0.8)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())  +
  geom_boxplot(alpha=0.6)
ponds_Shannon_18S $layers
ponds_Shannon_18S $layers <- ponds_Shannon_18S $layers[-1]
ponds_Shannon_18S <- ponds_Shannon_18S + scale_x_discrete(limits=summed_order) + scale_color_manual(values = region_palette)
ponds_Shannon_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("ponds_Shannon_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("ponds_Shannon_18S.pdf", width = 297, height = 210, units = c("mm"))

ponds_InvSimpson_18S <- plot_richness(ps_18S, x="new_pond", color ="region", shape = "water_body", measures=c("InvSimpson"),  title = "Inverse Simpson diversity index of 18S samples diversity grouped by sample sites")+
  geom_point(size=8, alpha=0.8)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())  +
  geom_boxplot(alpha=0.6)
ponds_InvSimpson_18S $layers
ponds_InvSimpson_18S $layers <- ponds_InvSimpson_18S $layers[-1]
ponds_InvSimpson_18S <- ponds_InvSimpson_18S + scale_x_discrete(limits=summed_order) + scale_color_manual(values = region_palette)
ponds_InvSimpson_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("ponds_InvSimpson_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("ponds_InvSimpson_18S.pdf", width = 297, height = 210, units = c("mm"))

#PLOT 4. Alpha diversity by regions
region_ACE_18S <- plot_richness(ps_18S, x="region", color ="region", measures=c("ACE"), title = "ACE diversity index of 18S samples grouped by region")+
  geom_jitter(size=8, alpha=0.8, aes(shape  = water_body), width = 0.2)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())
region_ACE_18S $layers
region_ACE_18S $layers <- region_ACE_18S $layers[-1]
region_ACE_18S <- region_ACE_18S + scale_color_manual(values = region_palette) + geom_boxplot(width=0.1, alpha=0.8)
region_ACE_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("region_ACE_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("region_ACE_18S.pdf", width = 297, height = 210, units = c("mm"))

region_Shannon_18S <- plot_richness(ps_18S, x="region", color ="region", measures=c("Shannon"), title = "Shannon diversity index of 18S samples grouped by region")+
  geom_jitter(size=8, alpha=0.8, aes(shape  = water_body), width = 0.2)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())
region_Shannon_18S $layers
region_Shannon_18S $layers <- region_Shannon_18S $layers[-1]
region_Shannon_18S <- region_Shannon_18S + scale_color_manual(values = region_palette) + geom_boxplot(width=0.1, alpha=0.8)
region_Shannon_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("region_Shannon_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("region_Shannon_18S.pdf", width = 297, height = 210, units = c("mm"))

region_InvSimpson_18S <- plot_richness(ps_18S, x="region", color ="region", measures=c("InvSimpson"), title = "Inverse Simpson diversity index of 18S samples grouped by region")+
  geom_jitter(size=8, alpha=0.8, aes(shape  = water_body), width = 0.2)+
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank())
region_InvSimpson_18S $layers
region_InvSimpson_18S $layers <- region_InvSimpson_18S $layers[-1]
region_InvSimpson_18S <- region_InvSimpson_18S + scale_color_manual(values = region_palette) + geom_boxplot(width=0.1, alpha=0.8)
region_InvSimpson_18S + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#use ggsave to save to your wd as a pdf and jpeg
ggsave("region_InvSimpson_18S.jpg", width = 297, height = 210, units = c("mm"))
ggsave("region_InvSimpson_18S.pdf", width = 297, height = 210, units = c("mm"))

#calculating richness and diversity####
#Null hypothesis = no difference in richness or evenness between sample regions.

#We will use the estimate_richness table we made earlier with the latitude and water temperature
richness_18S 

#We will add the region and new_ID metadata
richness_18S <- richness_18S %>% cbind(sample_data(ps_18S)[,"region"])%>% cbind(sample_data(ps_18S)[,"new_ID"])
richness_18S

#Save this output as a .txt file
write.table(richness_18S, file = "C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/dada2/final_results/richness_estimates_18S.txt")


#Check skew of data (normality)
hist(richness_18S$ACE)
hist(richness_18S$Shannon)
hist(richness_18S$InvSimpson)

#Test for normality (If p >0.05, data is normal)
shapiro.test(richness_18S$ACE) #p = 0.2174
shapiro.test(richness_18S$Shannon) #p = 3.541e-06
shapiro.test(richness_18S$InvSimpson) #p = 8.449e-09
#This shows us that the ACE data is normal but the evenness values for Shannon and Simpson are not normal

#As such we will use the non-parametric test for all diversity measures for continuity 

kruskal.test(richness_18S$ACE ~ richness_18S$region, data = richness_18S)
# data:  richness_18S$ACE by richness_18S$region
#Kruskal-Wallis chi-squared = 45.479, df = 6, p-value = 3.759e-08

ACE_18S_wilcoxon_summary <- pairwise.wilcox.test(richness_18S$ACE, richness_18S$region, p.adjust.method = "none")
ACE_18S_wilcoxon_summary
#                   1_kuujjuarapik 2_umiujaq
#2_umiujaq          0.65008        -        
#3_cambridge_bay    0.00014        0.00127  
#4_bylot_island     0.75830        0.93098  
#5_resolute         1.3e-05        0.00171  
#6_ellesmere_island 3.4e-06        0.00047  
#7_ice_shelves      3.4e-05        0.00355  
#                   3_cambridge_bay 4_bylot_island
#2_umiujaq          -               -             
#3_cambridge_bay    -               -             
#4_bylot_island     0.16529         -             
#5_resolute         0.01186         0.01732       
#6_ellesmere_island 1.1e-05         0.00433       
#7_ice_shelves      0.01252         0.03175       
#                   5_resolute 6_ellesmere_island
#2_umiujaq          -          -                 
#3_cambridge_bay    -          -                 
#4_bylot_island     -          -                 
#5_resolute         -          -                 
#6_ellesmere_island 0.06494    -                 
#7_ice_shelves      1.00000    0.09958    

kruskal.test(richness_18S$Shannon ~ richness_18S$region, data = richness_18S)
#data:  richness_18S$Shannon by richness_18S$region
#Kruskal-Wallis chi-squared = 22.98, df = 6, p-value = 0.000803

Shannon_18S_wilcoxon_summary <- pairwise.wilcox.test(richness_18S$Shannon, richness_18S$region, p.adjust.method = "none")
Shannon_18S_wilcoxon_summary
#data:  richness_18S$Shannon and richness_18S$region 

#                   1_kuujjuarapik 2_umiujaq
#2_umiujaq          0.33053        -        
#3_cambridge_bay    0.26963        0.09271  
#4_bylot_island     0.24509        0.48230  
#5_resolute         0.49394        0.27283  
#6_ellesmere_island 0.00100        0.00022  
#7_ice_shelves      0.05096        0.01295  
#                   3_cambridge_bay 4_bylot_island
#2_umiujaq          -               -             
#3_cambridge_bay    -               -             
#4_bylot_island     0.08629         -             
#5_resolute         0.80618         0.12554       
#6_ellesmere_island 0.00023         0.00433       
#7_ice_shelves      0.01905         0.00794       
#                   5_resolute 6_ellesmere_island
#2_umiujaq          -          -                 
#3_cambridge_bay    -          -                 
#4_bylot_island     -          -                 
#5_resolute         -          -                 
#6_ellesmere_island 0.00866    -                 
#7_ice_shelves      0.05195    0.08225  

kruskal.test(richness_18S$InvSimpson ~ richness_18S$region, data = richness_18S)
#data:  richness_18S$InvSimpson by richness_18S$region
#Kruskal-Wallis chi-squared = 20.044, df = 6, p-value = 0.00272

InvSimpson_18S_wilcoxon_summary <- pairwise.wilcox.test(richness_18S$InvSimpson, richness_18S$region, p.adjust.method = "none")
InvSimpson_18S_wilcoxon_summary

#data:  richness_18S$InvSimpson and richness_18S$region 

#                   1_kuujjuarapik 2_umiujaq
#2_umiujaq          0.11531        -        
#3_cambridge_bay    0.80225        0.11143  
#4_bylot_island     0.12902        0.41423  
#5_resolute         0.85975        0.52670  
#6_ellesmere_island 0.00510        0.00022  
#7_ice_shelves      0.38233        0.03747  
#                   3_cambridge_bay 4_bylot_island
#2_umiujaq          -               -             
#3_cambridge_bay    -               -             
#4_bylot_island     0.05505         -             
#5_resolute         0.33641         0.24675       
#6_ellesmere_island 0.00066         0.00866       
#7_ice_shelves      0.09922         0.01587       
#                   5_resolute 6_ellesmere_island
#2_umiujaq          -          -                 
#3_cambridge_bay    -          -                 
#4_bylot_island     -          -                 
#5_resolute         -          -                 
#6_ellesmere_island 0.00433    -                 
#7_ice_shelves      0.05195    0.12554 


