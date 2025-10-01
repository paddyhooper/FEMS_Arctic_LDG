#Making bar graphs for CoNet analysis
#Author: Patrick M. Hooper
#Created: 21.12.22

#Libraries
library(ggplot2); packageVersion("ggplot2")

#Network 1 - Prokaryotes and Protists and Fungi
net_1_proks <- read.table("~/16S_protist_fungi_proks.txt", header = TRUE)
net_1_euks <- read.table("~/16S_protist_fungi_euks.txt", header = TRUE)
net_1_proks
net_1_euks

#network 2 - 18S only
net_2_euks <- read.table("~/18S_only_euks.txt", header = TRUE)
net_2_euks

#Plotting Network 1 Graphs
#We should now be able to use GGPLOT to make our bar graphs, with the colour filled in with the col name column

p_net_1_proks <- ggplot(data = net_1_proks, aes(x=Degree, y=fct_reorder(Tax_Group, Degree), fill=I(col)))+geom_bar(stat="identity", width = 0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p_net_1_proks

#use ggsave to save to your wd as a pdf and jpeg
ggsave("net_1_proks_barplot.jpg", width = 297, height = 210, units = c("mm"))
ggsave("net_1_proks_barplot.pdf", width = 297, height = 210, units = c("mm"))

p_net_1_euks <- ggplot(data = net_1_euks, aes(x=Degrees, y=fct_reorder(Division, Degrees), fill=I(col)))+geom_bar(stat="identity", width = 0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p_net_1_euks

#use ggsave to save to your wd as a pdf and jpeg
ggsave("net_1_euks_barplot.jpg", width = 297, height = 210, units = c("mm"))
ggsave("net_1_euks_barplot.pdf", width = 297, height = 210, units = c("mm"))

#Plotting network 2 graphs
p_net_2_euks <- ggplot(data = net_2_euks, aes(x=Degrees, y=fct_reorder(Division, Degrees), fill=I(col)))+geom_bar(stat="identity", width = 0.6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p_net_2_euks

#use ggsave to save to your wd as a pdf and jpeg
ggsave("net_2_euks_barplot.jpg", width = 210, height = 297, units = c("mm"))
ggsave("net_2_euks_barplot.pdf", width = 210, height = 297, units = c("mm"))


