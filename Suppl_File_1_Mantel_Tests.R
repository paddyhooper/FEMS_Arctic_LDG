#Mantel Tests for 16S and 18S rRNA gene datasets
#Author: Patrick M. Hooper
#Date Created: 16/08/22

#REFERENCE:
#https://jkzorz.github.io/2019/07/08/mantel-test.html

#Mantel Tests are correlation tests that determine the correlation between two matrices (rather than two variables)
#In microbial ecology, the matrices are often distance matrices with corresponding positions
#They are a good way of comparing continous variables to our ASV counts

#To calculate the correlation, the matrix values of both matrices are unfolded into long column vectors, which are then used to determine correlation. https://sites.google.com/site/mb3gustame/hypothesis-tests/the-mantel-test

#To conduct a Mantel Test we need
#1. Species abundance dissimilarity matrix - created using a distance measure, i.e. Bray-Curtis dissimilarity
#2. Environmental parameter distance matrix - Generally created using Euclidean Distance (i.e. temeprature between samples)
#3. Geographic distance matrix - physical distances between sites

#With these we can determine if differences in community composition between samples are correlated or co-vary with the differences in temperature between samples or the physical distance between samples, i.e. is the environment selecting or is there a geographical limitation on diversity across the gradient?

set.seed(2305)

#STEP 1. LOAD PACKAGES#### 
library(vegan); packageVersion("vegan")
library(microbiome); packageVersion("microbiome")
library(geosphere); packageVersion("geosphere")

#STEP 2. LOAD YOUR AVERAGED PHYLOSEQ OBJECTS FROM THE ORDINATION SCRIPT####
ps_16S_ave <- readRDS(file = "~/ps_16S_ave")
ps_func_18S_ave <- readRDS(file = "~/ps_func_18S_ave")

#STEP 3. RELATIVE ABUNDANCE TRANSFORM THE TWO PHYLOSEQ OBJECTS TO MATCH THE ORDINATION STEPS####
ps_16S_ave_RA = transform_sample_counts(ps_16S_ave, function(x) x / sum(x))
ps_func_18S_ave_RA  = transform_sample_counts(ps_func_18S_ave, function(x) x / sum(x))

#STEP 4. LOAD YOUR COUNT DATA FOR 16s####
count_16S_ra <- abundances(ps_16S_ave_RA)
count_16s_ra <- as.data.frame(count_16S_ra)
#This needs to be transposed so the samples are rows
count_16S_ra <- t(count_16S_ra)
rownames(count_16S_ra) # check your samples are rows

#STEP 5. MAKE A TABLE OF COUNT DATA FOR 18s####
count_18S_ra <- abundances(ps_func_18S_ave_RA)
count_18S_ra <- as.data.frame(count_18S_ra)
#This needs to be transposed so the samples are rows
count_18S_ra <- t(count_18S_ra)
rownames(count_18S_ra) # check your samples are rows

#STEP 6. LOAD YOUR 16S AND 18S METADATA TABLE####
#Load the 16S average metadata table
meta_average_16S <- read.table(file = "~/GLOM_V4_16S_METADATA_FORMATTED_AVERAGE_UPDATE_17_8_22.txt")
nrow(meta_average_16S)# 32 samples
colnames(meta_average_16S)
#Subset the metadata into smaller dataframes with the necessary information

#Load the 18S average metadata table
meta_average_18S <- read.table(file = "~/GLOM_V9_18S_METADATA_FORMATTED_UPDATE_AVERAGE_17_8_22.txt")
nrow(meta_average_18S)# 32 samples
colnames(meta_average_18S)

#You're all set!

#MANTEL TEST - 16S####
#Subset the metadata into smaller dataframes with the necessary information
#environmental vector: water temp, pH, conductivity and average annual temperature
env_16S <- as.data.frame(meta_average_16S)
env_16S <- env_16S[,13:17]
env_16S

#scale data - environmental data needs to be scaled prior to creating a data matrix as the environmental variables were all measured using different metrics
scale.env_16S <- scale(env_16S, center = TRUE, scale = TRUE)
scale.env_16S #take a look

#Long lat data
geo_16S <- as.data.frame(meta_average_16S)
geo_16S <- geo_16S[,9:10]
geo_16S <- geo_16S[,2:1] #swap long and lat around to make geosphere happy
geo_16S

#Now we have to convert these subsets into distance matrices
#abundance data frame - bray curtis dissimilarity
dist.abund_16S <- vegdist(count_16S_ra, method = "bray")

#create distance matrix of scaled data of water pH, conductivity, temp, and air temperature to see influence of environmental factors
dist.env_16S <- dist(scale.env_16S, method = "euclidean")

#We then look at each of these variables one by one
#water temp vector - euclidean distance
dist.temp_16S <- dist(env_16S$water_temp, method = "euclidean")

#air temp vector - euclidean distance
dist.air.temp_16S <- dist(env_16S$ave_ann_temp, method = "euclidean")

#sample year temp vector - euclidean distance
dist.sample.year.temp_16S <- dist(env_16S$sample_year_temp, method = "euclidean")

#water pH - euclidean distance
dist.pH_16S <- dist(env_16S$water_pH, method = "euclidean")

#water conductivity - euclidean distance
dist.cond_16S <- dist(env_16S$water_conductivity, method = "euclidean")

#geographic data frame - haversine distance 
#The Haversine (or great circle) distance is the angular distance between two points on the surface of a sphere
d.geo_16S <- distm(geo_16S, fun = distHaversine)
dist.geo_16S <- as.dist(d.geo_16S)

#Now we run the mantel command
#The mantel command requires the user to specify certain parameters:
#1. distance matrices (i.e. dist.abund and dist.temp)
#2. correlation method. I use Spearman to make the test "non-parametric". Learn more about correlation methods here permutations. Mantel tests determine significance by permuting (randomizing) one matrix X number of times and observing the expected distribution of the statistic. I tend to pick a larger permutation number, but if computing power is low, feel free to decrease 
#3. na.rm. An optional addition to the command that tells R to delete rows in which there are missing values.

#what does permutations mean? 
#Mantel tests determine significance by permuting (randomizing) one matrix X number of times and observing the expected distribution of the statistic. I tend to pick a larger permutation number, but if computing power is low, feel free to decrease

#16s - abundance vs environmental 
abund_env_16S <- mantel(dist.abund_16S, dist.env_16S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_env_16S

#Mantel statistic r: 0.162 
#Significance: 0.0124 

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0793 0.1084 0.1358 0.1676 
#Permutation: free
#Number of permutations: 9999

#THERE IS A SIGNIFICANT CORRELATION BETWEEN ALL ENVIRONMENTAL VARIABLES (AIR TEMP, WATER TEMP, WATER PH, WATER CONDUCTIVITY AND bray-curtis dissimilarity )

#16s - abundance vs water temp
abund_water_temp_16S <- mantel(dist.abund_16S, dist.temp_16S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_water_temp_16S

#Mantel statistic r: -0.04954 
#Significance: 0.7132 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.104 0.142 0.176 0.217 
#Permutation: free
#Number of permutations: 9999
#THERE IS NO SIGNIFICANT CORRELATION BETWEEN  bray-curtis dissimilarity AND WATER TEMPERATURE

#16s - abundance vs annual air temp
abund_ann_air_temp_16S <- mantel(dist.abund_16S, dist.air.temp_16S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_ann_air_temp_16S

#Mantel statistic r: 0.2141 
#Significance: 0.0017 

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0685 0.0962 0.1237 0.1576 
#Permutation: free
#Number of permutations: 9999

#THERE IS SIGNIFICANT CORRELATION BETWEEN AIR TEMPERATURE AND  bray-curtis dissimilarity COUNT IN 16S

#16s - abundance vs sample year temp
abund_sample_year_temp_16S <- mantel(dist.abund_16S, dist.sample.year.temp_16S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_sample_year_temp_16S 

#Mantel statistic r: 0.2321 
#Significance: 0.0013 

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0678 0.0917 0.1178 0.1548 
#Permutation: free
#Number of permutations: 9999

#There is a significant correlation between sample year temp and 16S distance

#16s - abundance vs pH
abund_pH_16S <- mantel(dist.abund_16S, dist.pH_16S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_pH_16S

#Mantel statistic r: 0.09097 
#Significance: 0.0695 

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0756 0.1042 0.1309 0.1647 
#Permutation: free
#Number of permutations: 9999

#THERE IS NO SIGNIFICANT CORRELATION BETWEEN PH AND  bray-curtis dissimilarity  COUNT IN 16S

#16s - abundance vs conductivity
abund_cond_16S <- mantel(dist.abund_16S, dist.cond_16S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_cond_16S

#Mantel statistic r: 0.002144 
#Significance: 0.4446 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.110 0.148 0.183 0.221 
#Permutation: free
#Number of permutations: 9999

#THERE IS NO SIGNIFICANT CORRELATION BETWEEN CONDUCTIVITY AND  bray-curtis dissimilarity  IN 16S

#16s - abundance vs geographic (vs latitude and longitude)
abund_geo_16S <- mantel(dist.abund_16S, dist.geo_16S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo_16S

#Mantel statistic r: 0.2083 
#Significance: 0.002

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0671 0.0932 0.1204 0.1482
#Permutation: free
#Number of permutations: 9999

#THERE IS A STRONG SIGNIFICANT CORRELATION BETWEEN DISTANCE AND  bray-curtis dissimilarity  ON THE INDIVIDUAL SAMPLES

#Mantel test on 18S relative abundance data####
#Firstly, we need a table of abundance data
#Utilise the average phyloseq object we made in our ordination plot

#environmental vector: water temp, pH, conductivity and average annual temperature
env_18S <- as.data.frame(meta_average_18S)
env_18S <- env_18S[,13:17]
env_18S

#scale data - environmental data needs to be scaled prior to creating a data matrix as the environmental variables were all measured using different metrics
scale.env_18S <- scale(env_18S, center = TRUE, scale = TRUE)

#Long lat data
geo_18S <- as.data.frame(meta_average_18S)
geo_18S <- geo_18S[,9:10]
geo_18S <- geo_18S[,2:1] #swap long and lat around to make geosphere happy
geo_18S

#Now we have to convert these subsets into distance matrices
#abundance data frame - bray curtis dissimilarity
dist.abund_18S <- vegdist(count_18S_ra, method = "bray")

#create distance matrix of scaled data of water pH, conductivity, temp, and air temperature to see influence of environmental factors
dist.env_18S <- dist(scale.env_18S, method = "euclidean")

#We then look at each of these variables one by one
#water temp vector - euclidean distance
dist.temp_18S <- dist(env_18S$water_temp, method = "euclidean")

#air temp vector - euclidean distance
dist.air.temp_18S <- dist(env_18S$ave_ann_temp, method = "euclidean")

#sample year temp vector - euclidean distance
dist.sample.year.temp_18S <- dist(env_18S$sample_year_temp, method = "euclidean")

#water pH - euclidean distance
dist.pH_18S <- dist(env_18S$water_pH, method = "euclidean")

#water conductivity - euclidean distance
dist.cond_18S <- dist(env_18S$water_conductivity, method = "euclidean")

#geographic data frame - haversine distance 
d.geo_18S <- distm(geo_18S, fun = distHaversine)
dist.geo_18S <- as.dist(d.geo_18S)

#Now we run the mantel command
#The mantel command requires the user to specify certain parameters:
#1. distance matrices (i.e. dist.abund and dist.temp)
#2. correlation method. I use Spearman to make the test "non-parametric". Learn more about correlation methods here permutations. Mantel tests determine significance by permuting (randomizing) one matrix X number of times and observing the expected distribution of the statistic. I tend to pick a larger permutation number, but if computing power is low, feel free to decrease 
#3. na.rm. An optional addition to the command that tells R to delete rows in which there are missing values.

#what does permutations mean? 
#Mantel tests determine significance by permuting (randomizing) one matrix X number of times and observing the expected distribution of the statistic. I tend to pick a larger permutation number, but if computing power is low, feel free to decrease

#abundance vs environmental 
abund_env_18S <- mantel(dist.abund_18S, dist.env_18S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_env_18S

#Mantel statistic r: 0.07279 
#Significance: 0.1654 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.101 0.134 0.161 0.195 
#Permutation: free
#Number of permutations: 9999

#No significant association between environment and  bray-curtis dissimilarity 

#abundance vs water temp
abund_water_temp_18S <- mantel(dist.abund_18S, dist.temp_18S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_water_temp_18S

#Mantel statistic r: -0.1535 
#Significance: 0.9522 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.136 0.183 0.226 0.275 
#Permutation: free
#Number of permutations: 9999

#No significant association between  bray-curtis dissimilarity  and water temperature

#abundance vs air temp
abund_air_temp_18S <- mantel(dist.abund_18S, dist.air.temp_18S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_air_temp_18S

#Mantel statistic r: 0.117 
#Significance: 0.031

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0749 0.1004 0.1233 0.1502 
#Permutation: free
#Number of permutations: 9999

#Significant association between air temperature and  bray-curtis dissimilarity 

#abundance vs sample year temp
abund_sample_year_temp_18S <- mantel(dist.abund_18S, dist.sample.year.temp_18S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_sample_year_temp_18S 

#Mantel statistic r: 0.1525 
#Significance: 0.0102 

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0759 0.1030 0.1264 0.1528 
#Permutation: free
#Number of permutations: 9999

#There is a significant correlation between sample year temperature and dissimilarity in the 18S table

#abundance vs pH
abund_pH_18S <- mantel(dist.abund_18S, dist.pH_18S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_pH_18S

#Mantel statistic r: 0.121 
#Significance: 0.0462

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0877 0.1169 0.1449 0.1707 
#Permutation: free
#Number of permutations: 9999

#Significant association between pH  bray-curtis dissimilarity 

#abundance vs conductivity
abund_cond_18S <- mantel(dist.abund_18S, dist.cond_18S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_cond_18S

#Mantel statistic r: -0.09993 
#Significance: 0.8416 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.139 0.181 0.223 0.267 
#Permutation: free
#Number of permutations: 9999

#No significant association between conductivity and bray-curtis dissimilarity 

#abundance vs geographic (vs latitude and longitude)
abund_geo_18S <- mantel(dist.abund_18S, dist.geo_18S, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo_18S

#Mantel statistic r: 0.1102 
#Significance: 0.0396 

#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0766 0.1020 0.1250 0.1588 
#Permutation: free
#Number of permutations: 9999

#Significant association between geographic distance and bray-curtis dissimilarity 

