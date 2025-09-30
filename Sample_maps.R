#Mapping sampling regions in Canada
#Author: Patrick M Hooper
#Date created: 22/11/22

#Adapted from - Plotting in ggplot2
#https://r-spatial.org/r/2018/10/25/ggplot2-sf.html

#This script produces a map of the sampling locations for the latitudinal diversity gradient study
#Four maps will be produced:
  #1. Total sampling region
  #2. KJ and UM
  #3. CB and BI and RI
  #4. High Arctic

#setwd
setwd("C:/Users/pmh36/OneDrive - Natural History Museum/R/R_data/map")

#Install packages
install.packages(c("rnaturalearth", "rnaturalearthdata", "ggspatial"))

#Load packages
library("tidyverse")
library("ggplot2")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library("ggrepel")

#ggplot
theme_set(theme_bw())

#1. Total sample region
#create base map from world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#add latitudinal limits for the total sample region
map <- ggplot(data = world) +
  geom_sf() + annotation_scale(location = "bl", width_hint = 0.5) + coord_sf(xlim = c(-115, -60), ylim = c(50, 85), expand = FALSE)
map

#add the style points
format_map <- map + theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) + xlab("Longitude") + ylab("Latitude") + ggtitle("Sampling sites across the Canadian Arctic")
format_map

#Load a .csv file of the sampling points consisting of latitude and longitude
#NOTE: must be decimal degrees in R
samps = read.csv(file = "sample_regions.csv", fileEncoding="UTF-8-BOM")
samps

#Add sampling points to map
format_map + geom_point(data = samps, aes(x = lon, y = lat), size = 5, shape = samps$shape, col = alpha(samps$col, 0.8))

#Add your sample region names to your map
locs = read.csv(file = "sampling_regions.csv", fileEncoding="UTF-8-BOM")

format_map + geom_point(data = samps, aes(x = lon, y = lat), size = 5, shape = samps$shape, col = alpha(samps$col, 0.8)) + geom_text_repel(data = locs, aes(x = lon, y = lat, label = id))

#Save your map
ggsave("total_sample_region.pdf", width = 210, height = 297, units = c("mm"))
ggsave("total_sample_region.png", width = 210, height = 297, units = c("mm"))

#Now we have the basic formula, we just need to change the latitude and longitude scale to create subsets for each sample region

#2. Taiga Samples
#add latitudinal limits for KJ and UM samples

#Subset only Kuujjuarapik and Umiujaq labels from locs dataframe
taiga_locs <- locs[1:3,]
taiga_locs

ggplot(data = world) +
  geom_sf() + annotation_scale(location = "bl", width_hint = 0.5) + coord_sf(xlim = c(-70, -90), ylim = c(50, 60), expand = FALSE)+ theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) + xlab("Longitude") + ylab("Latitude") + ggtitle("Map of sample locations")+ geom_point(data = samps, aes(x = lon, y = lat), size = 5, shape = samps$shape, col = alpha(samps$col, 0.8)) + geom_text_repel(data = taiga_locs, aes(x = lon, y = lat, label = id))


