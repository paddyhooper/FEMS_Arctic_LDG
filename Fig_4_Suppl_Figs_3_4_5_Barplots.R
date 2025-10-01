#16S and 18S Barplots
#Author: Patrick M Hooper
#Date Created: 23/12/22

#This script was used to make the taxonomic and functional plots for the 16S and 18S data
#Functional annotation was added to a taxonomy table based on literature review. This information was added to the phyloseq object created in the previous script '16S_18S_ordination.R'
#This script uses the taxonomic data and functional data from the saved phyloseq objects for the 16S and 18S data

#Set up####
#1. Load Packages
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(Biostrings); packageVersion("Biostrings")
library(microbiome); packageVersion("microbiome")
library(ampvis2); packageVersion("ampvis2")
library(forcats); packageVersion("forcats")
library(chroma); packageVersion("chroma")
library(viridis); packageVersion("viridis")

#2. Set ggplot2 theme
theme_set(theme_minimal()) 

#Set standardized plot font sizes
font_size <- theme(axis.text = element_text(size = 16)) + theme(axis.title = element_text(size = 16)) + theme(plot.title = element_text(size = 20)) + theme(legend.text = element_text(size = 16)) + theme(legend.title = element_text(size = 16)) + theme(plot.subtitle = element_text(size = 18))

#Load your data####
#4. Load your phyloseq objects from saveRDS
#These files were created in script '16S_18S_ordination.R'
ps_16S_ave <- readRDS(file = "~/ps_16S_ave")
ps_16S_ave #15250 taxa
ps_func_18S_ave <- readRDS(file = "~/ps_func_18S_ave")
ps_func_18S_ave #9354 taxa

#5. Relative abundance transformation of the count data
ps_16S_ave_RA = transform_sample_counts(ps_16S_ave, function(x) x / sum(x))
ps_func_18S_ave_RA  = transform_sample_counts(ps_func_18S_ave, function(x) x / sum(x))

#6. Create a function with your samples ordered by latitude from lowest to highest
average_order = c("KJ1","KJ2","KJ3","KJ4","KJ5","KJ6","KJ7","KJ8","UM1","UM2","UM3","UM4","UM5","UM6","UM7","UM8","CB1","CB2","CB3","CB4","CB5","CB6","CB7","CB8","BY1","BY2","RE1","RE2","WH","AP","MKIS","WHIS")

#7. Set your figure colour palette
#We can then extend this for the different plots using the chroma package
figure_palette <- c("#228833", "#CCBB44", "#EE6677", "#AA3377", "#66CCEE", "#4477AA", "#BBBBBB")
show_col(figure_palette)

#Plotting taxonomic data####
#This section we will create bubble plots of the 16S and 18S data

#16S taxonomic plot####
#8. Prepare you 16S data

#Make a bubble plot of the 16S data at the phylum level with the low abundance ASVS grouped
#This requires us collating the low abundance ASVs (determined here at <0.1% average abundance)
phylum_glom <- tax_glom(ps_16S_ave_RA, taxrank = 'Phylum')
phylum_glom #49 phyla
phylum_glom_melt <- psmelt(phylum_glom) #create a dataframe you can use in ggplot2
phylum_glom_melt$Phylum <- as.character(phylum_glom_melt$Phylum) #convert the Phylum column to a character

#Group the low abundance data (<1%)
phylum_glom_melt$Phylum[phylum_glom_melt$Abundance < 0.01] <- "<1% total abundance"

#How many phyla are you left with?
unique(phylum_glom_melt$Phylum) #This leaves us with 18 phyla and 1 grouped abundance

#Make an extended colour palette with this many colours from your figure palette
prok_palette <- interp_colors(20, colors = figure_palette)
prok_palette

#9. Plot a bubble plot
bubble_plot_16S_phyla <- ggplot(phylum_glom_melt, aes(x=new_pond, y = fct_reorder(Phylum, Abundance))) + geom_point(aes(size = Abundance, colour = fct_rev(fct_reorder(Phylum, Abundance))), alpha = 0.8)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,25), breaks = c(0.10, 0.25, 0.50, 0.75, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_manual(values = prok_palette, guide = "none")+ scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())
bubble_plot_16S_phyla + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + expand_limits(y = 20)

#Use ggsave to save a pdf and jpeg to the working directory in A4 ratio
ggsave("bubble_plot_16S_phyla.jpg", width = 297, height = 210, units = c("mm"))
ggsave("bubble_plot_16S_phyla.pdf", width = 297, height = 210, units = c("mm"))

#18S taxonomic plot####
#10. Prepare your 18S data
division_glom <- tax_glom(ps_func_18S_ave_RA, taxrank = 'Division')
division_glom #29 phyla
division_glom_melt <- psmelt(division_glom)
division_glom_melt$Division <- as.character(division_glom_melt$Division)

#Group the low abundance data (<0.01%)
division_glom_melt$Division[division_glom_melt$Abundance < 0.001] <- "<0.1% total abundance"

#How many divisions are you left with?
unique(division_glom_melt$Supergroup)

#Make an extended colour palette with this many colours from your figure palette
euk_palette <- interp_colors(9, colors = figure_palette)
euk_palette

#11. Plot a bubble plot
bubble_plot_18S_division <- ggplot(division_glom_melt, aes(x=new_pond, y = fct_reorder(Division, Abundance))) + geom_point(aes(size = Abundance, colour = fct_rev(fct_reorder(Supergroup, Abundance))), alpha = 0.8,)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,15), breaks = c(0.10, 0.25, 0.50, 0.75, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_manual(values = euk_palette)+ scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())
bubble_plot_18S_division+ font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + expand_limits(y = 20) 

#Save your plot
ggsave("bubble_plot_18S_division.jpg", width = 297, height = 210, units = c("mm"))
ggsave("bubble_plot_18S_division.pdf", width = 297, height = 210, units = c("mm"))

                                              #Fig. 4. Cyanobacteria Plot###
#cyanobacteria <- subset_taxa(ps_16S_ave_RA, Phylum == "Cyanobacteria")
#cyanobacteria_glom <- tax_glom(cyanobacteria, taxrank = 'Family')
#cyanobacteria_glom #15 orders
#cyanobacteria_glom_melt <- psmelt(cyanobacteria_glom) #create a dataframe you can use in ggplot2
#cyanobacteria_glom_melt$Family <- as.character(cyanobacteria_glom_melt$Family) #convert the Phylum column to a character

#Group the low abundance data (<1%)
#cyanobacteria_glom_melt$Family[cyanobacteria_glom_melt$Abundance < 0.0001] <- "<0.01% total abundance"

#How many phyla are you left with?
#unique(cyanobacteria_glom_melt$Family) #This leaves us with 15 classes and 1 grouped abundance

#Make an extended colour palette with this many colours from your figure palette
#cyano_palette <- interp_colors(16, colors = figure_palette)
#cyano_palette

#9. Plot a bubble plot
#bubble_plot_cyano <- ggplot(cyanobacteria_glom_melt, aes(x=new_pond, y = fct_rev(fct_infreq(Family)))) + geom_point(aes(size = Abundance, colour = fct_rev(fct_reorder(Order, Abundance))), alpha = 0.8)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,25), breaks = c(0.10, 0.25, 0.50, 0.75, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_manual(values = cyano_palette)+ scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())
#bubble_plot_cyano+ font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + expand_limits(y = 2)

#Use ggsave to save a pdf and jpeg to the working directory in A4 ratio
#ggsave("bubble_plot_cyano.jpg", width = 297, height = 210, units = c("mm"))
#ggsave("bubble_plot_cyano.pdf", width = 297, height = 210, units = c("mm"))
                                            
#Fig. 4: Proteobacteria Plot####
#proteobacteria <- subset_taxa(ps_16S_ave_RA, Phylum == "Proteobacteria")
#proteobacteria_glom <- tax_glom(proteobacteria, taxrank = 'Order')
#proteobacteria_glom #15 orders
#proteobacteria_glom_melt <- psmelt(proteobacteria_glom) #create a dataframe you can use in ggplot2
#proteobacteria_glom_melt$Order <- as.character(proteobacteria_glom_melt$Order) #convert the Phylum column to a character

#Group the low abundance data (<1%)
#proteobacteria_glom_melt$Order[proteobacteria_glom_melt$Abundance < 0.0001] <- "<0.1% total abundance"

#How many phyla are you left with?
#unique(proteobacteria_glom_melt$Order) #This leaves us with 15 classes and 1 grouped abundance

#Make an extended colour palette with this many colours from your figure palette
#proteo_palette <- interp_colors(55, colors = figure_palette)
#proteo_palette

#9. Plot a bubble plot
#bubble_plot_proteo <- ggplot(proteobacteria_glom_melt, aes(x=new_pond, y = fct_rev(fct_infreq(Order)))) + geom_point(aes(size = Abundance, colour = fct_rev(fct_reorder(Order, Abundance))), alpha = 0.8)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,25), breaks = c(0.10, 0.25, 0.50, 0.75, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_manual(values = proteo_palette)+ scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())
#bubble_plot_proteo+ font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + expand_limits(y = 2)

#Use ggsave to save a pdf and jpeg to the working directory in A4 ratio
#ggsave("bubble_plot_proteo.jpg", width = 297, height = 210, units = c("mm"))
#ggsave("bubble_plot_proteo.pdf", width = 297, height = 210, units = c("mm"))

#SUPPL FIG. 3: Metazoa####
metazoa <- subset_taxa(ps_func_18S_ave_RA, Division == "Metazoa")
metazoa_glom <- tax_glom(metazoa, taxrank = 'Order')
metazoa_glom #24 orders
metazoa_glom_melt <- psmelt(metazoa_glom) #create a dataframe you can use in ggplot2
metazoa_glom_melt$Order <- as.character(metazoa_glom_melt$Order) #convert the Phylum column to a character

#Group the low abundance data (<1%)
metazoa_glom_melt$Order[metazoa_glom_melt$Abundance < 0.001] <- "<0.1% total abundance"

#How many phyla are you left with?
unique(metazoa_glom_melt$Order) #This leaves us with 15 orders and 1 grouped abundance

#Make an extended colour palette with this many colours from your figure palette
met_palette <- interp_colors(15, colors = figure_palette)
met_palette

#9. Plot a bubble plot
bubble_plot_metazoa <- ggplot(metazoa_glom_melt, aes(x=new_pond, y = fct_rev(fct_infreq(Order)))) + geom_point(aes(size = Abundance, colour = fct_rev(fct_reorder(Class, Abundance))), alpha = 0.8)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,25), breaks = c(0.10, 0.25, 0.50, 0.75, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_manual(values = met_palette)+ scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())
bubble_plot_metazoa + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + expand_limits(y = 2)

#Use ggsave to save a pdf and jpeg to the working directory in A4 ratio
ggsave("bubble_plot_metazoa.jpg", width = 297, height = 210, units = c("mm"))
ggsave("bubble_plot_metazoa.pdf", width = 297, height = 210, units = c("mm"))

#SUPPL FIG. 4: Fungi
fungi <- subset_taxa(ps_func_18S_ave_RA, Division == "Fungi")
fungi_glom <- tax_glom(fungi, taxrank = 'Order')
fungi_glom #19 orders
fungi_glom_melt <- psmelt(fungi_glom) #create a dataframe you can use in ggplot2
fungi_glom_melt$Order <- as.character(fungi_glom_melt$Order) #convert the Phylum column to a character

#Group the low abundance data (<1%)
fungi_glom_melt$Order[fungi_glom_melt$Abundance < 0.0001] <- "<0.01% total abundance"

#How many phyla are you left with?
unique(fungi_glom_melt$Order) #This leaves us with 15 classes and 1 grouped abundance

#Make an extended colour palette with this many colours from your figure palette
fungi_palette <- interp_colors(13, colors = figure_palette)
fungi_palette

#9. Plot a bubble plot
bubble_plot_fungi <- ggplot(fungi_glom_melt, aes(x=new_pond, y = fct_rev(fct_infreq(Order)))) + geom_point(aes(size = Abundance, colour = fct_rev(fct_reorder(Class, Abundance))), alpha = 0.8)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,25), breaks = c(0.10, 0.25, 0.50, 0.75, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_colour_manual(values = fungi_palette)+ scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())
bubble_plot_fungi+ font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1)) + expand_limits(y = 2)

#Use ggsave to save a pdf and jpeg to the working directory in A4 ratio
ggsave("bubble_plot_fungi.jpg", width = 297, height = 210, units = c("mm"))
ggsave("bubble_plot_fungi.pdf", width = 297, height = 210, units = c("mm"))

                                              
#SUPPL FIG. 5: Protist functional annotation####
#This step first requires filtering of all none protist groups from our 18S phyloseq object
ps_protist = subset_taxa(ps_func_18S_ave_RA, !Division %in% c("Fungi","Rhodophyta", "Metazoa") & !Class  %in% c("Embryophyceae"))
ps_protist #6296 taxa
protist_melt <- psmelt(ps_protist)

#Make colour palette
unique(protist_melt$Nutrition)

#Let's make a manual nutrition palette
#We need 11 colours, so let's extend our figure palette and find out the hex values
pal_11 <- (interp_colors(11, colors = figure_palette))
show_col(pal_11)
pal_11

function_palette <-
  c(
    "phototrophic" = "#228833",
    "plant_parasite" = "#91A73D",
    "bacterivore" = "#D5AC51",
    "omnivore" = "#E97A6F",
    "eukaryvore" = "#D35277",
    "mixotrophic" = "#AA3377",
    "osmotrophic" = "#9895BC",
    "animal_parasite" = "#5FBAE0",
    "pathogen" = "#4B87B7",
    "detritivore" = "#7991B1",
    "NA" = "#BBBBBB"
  ) 
show_col(function_palette)

#Make a bubble plot of the bacterivores####
bacterivores = subset_taxa(ps_protist, Nutrition %in% c("bacterivore")) 
bacterivores #1289 ASVs
class.sum = tapply(taxa_sums(bacterivores), tax_table(bacterivores)[, "Class"], sum, na.rm=TRUE)
top15_classes = names(sort(class.sum, TRUE))[1:15]
top15_classes
bacterivores = prune_taxa((tax_table(bacterivores)[, "Class"] %in% top15_classes), bacterivores)
bacterivores <- tax_glom(bacterivores, taxrank = "Class")
bacterivores_melt <- psmelt(bacterivores)

#Plot a bubble plot of the bacterivores, colour coordinated
bacterivores_plot <- ggplot(bacterivores_melt, aes(x=new_pond, y = fct_reorder(Class, Abundance))) + geom_point(aes(size = Abundance, colour = Kingdom), alpha = 0.8,)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,25), breaks = c(0.10, 0.25, 0.50, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())+scale_colour_manual(values = "#D5AC51", guide = "none")

bacterivores_plot + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#Save your plot
ggsave("bacterivores.jpg", width = 315, height = 297, units = c("mm"))
ggsave("bacterivores.pdf", width = 315, height = 297, units = c("mm"))

#Make a plot of the omnivores####
omnivores = subset_taxa(ps_protist, Nutrition %in% c("omnivore")) 
omnivores #471 ASVs
class.sum = tapply(taxa_sums(omnivores), tax_table(omnivores)[, "Class"], sum, na.rm=TRUE)
top15_classes = names(sort(class.sum, TRUE))[1:15]
top15_classes
omnivores = prune_taxa((tax_table(omnivores)[, "Class"] %in% top15_classes), omnivores)
omnivores <- tax_glom(omnivores, taxrank = "Class")
omnivores_melt <- psmelt(omnivores)

#Plot a bubble plot of the bacterivores, colour coordinated
omnivores_plot <- ggplot(omnivores_melt, aes(x=new_pond, y = fct_reorder(Class, Abundance))) + geom_point(aes(size = Abundance, colour = Kingdom), alpha = 0.8,)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,25), breaks = c(0.10, 0.25, 0.50, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())+scale_colour_manual(values = "#E97A6F", guide = "none")
omnivores_plot + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#Save your plot
ggsave("omnivores.jpg", width = 315, height = 297, units = c("mm"))
ggsave("omnivores.pdf", width = 315, height = 297, units = c("mm"))

#Make a plot of the eukaryvores
eukaryvores = subset_taxa(ps_protist, Nutrition %in% c("eukaryvore")) 
eukaryvores #420 ASVs
order.sum = tapply(taxa_sums(eukaryvores), tax_table(eukaryvores)[, "Order"], sum, na.rm=TRUE)
top15_order = names(sort(order.sum, TRUE))[1:15]
top15_order
eukaryvores = prune_taxa((tax_table(eukaryvores)[, "Order"] %in% top15_order), eukaryvores)
eukaryvores <- tax_glom(eukaryvores, taxrank = "Order")
eukaryvores_melt <- psmelt(eukaryvores)

#Plot a bubble plot of the bacterivores, colour coordinated
eukaryvores_plot <- ggplot(eukaryvores_melt, aes(x=new_pond, y = fct_reorder(Order, Abundance))) + geom_point(aes(size = Abundance, colour = Kingdom), alpha = 0.8,)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,25), breaks = c(0.10, 0.25, 0.50, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())+scale_colour_manual(values = "#D35277", guide = "none")
eukaryvores_plot + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#Save your plot
ggsave("eukaryvores.jpg", width = 315, height = 297, units = c("mm"))
ggsave("eukaryvores.pdf", width = 315, height = 297, units = c("mm"))

#Plot mixotrophic groups 
mixotrophic = subset_taxa(ps_protist, Nutrition %in% c("mixotrophic")) 
mixotrophic #315 ASVs
order.sum = tapply(taxa_sums(mixotrophic), tax_table(mixotrophic)[, "Order"], sum, na.rm=TRUE)
top15_order = names(sort(order.sum, TRUE))[1:15]
top15_order
mixotrophic = prune_taxa((tax_table(mixotrophic)[, "Order"] %in% top15_order), mixotrophic)
mixotrophic <- tax_glom(mixotrophic, taxrank = "Order")
mixotrophic_melt <- psmelt(mixotrophic)

#Plot a bubble plot of the bacterivores, colour coordinated
mixotrophic_plot <- ggplot(mixotrophic_melt, aes(x=new_pond, y = fct_reorder(Order, Abundance))) + geom_point(aes(size = Abundance, colour = Kingdom), alpha = 0.8,)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,25), breaks = c(0.10, 0.25, 0.50, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())+scale_colour_manual(values = "#AA3377", guide = "none")
mixotrophic_plot + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#Save your plot
ggsave("mixotrophic.jpg", width = 315, height = 297, units = c("mm"))
ggsave("mixotrophic.pdf", width = 315, height = 297, units = c("mm"))

#Plot animal parasites 
parasites = subset_taxa(ps_protist, Nutrition %in% c("animal_parasite"))
parasites #275 taxa
order.sum = tapply(taxa_sums(parasites), tax_table(parasites)[, "Order"], sum, na.rm=TRUE)
top15_order = names(sort(order.sum, TRUE))[1:15]
top15_order
parasites = prune_taxa((tax_table(parasites)[, "Order"] %in% top15_order), parasites)
parasites <- tax_glom(parasites, taxrank = "Order")
parasites_melt <- psmelt(parasites)

#Plot a bubble plot of the bacterivores, colour coordinated
parasites_plot <- ggplot(parasites_melt, aes(x=new_pond, y = fct_reorder(Order, Abundance))) + geom_point(aes(size = Abundance, colour = Kingdom), alpha = 0.8,)+ scale_size_continuous(limits = c(0.00001, 1), range = c(1,25), breaks = c(0.10, 0.25, 0.50, 1))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_x_discrete(limits=average_order) + labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") + theme(strip.background = element_blank())+scale_colour_manual(values = "#5FBAE0", guide = "none")
parasites_plot + font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))

#Save your plot
ggsave("parasite_plot.jpg", width = 315, height = 297, units = c("mm"))
ggsave("parasite_plot.pdf", width = 315, height = 297, units = c("mm"))

#I didn't make a plot for osmotrophic, plant_parasite, pathogens, or detritivores as they were all represented by one division, will be discussed instead. 

#### END OF SCRIPT ####

