###These are the next steps in calculating genetic diversity 
##First, we tested linkage disequilibrium, null alleles, and Hardy Weinberg 
#equilibrium. Next, a data frame including expected heterozygosity, 
#allelic richness, number of alleles, mean longtitude and latitude 
#for wild populations, and individual numbers. 
#This table is included in full in the supplemental text of this manuscript.
#When files are referred to as "clean" that means individuals 
#that are clones and indviduals with too much missing data have been removed. 
#When files and objects are titled "red" that means they have been reduced
#for relatedness (25% or more related individuals are reduced to one individual
#per phenotype)

#########################
#        Libraries      #
#########################

library(adegenet)
library(poppr)
library(hierfstat)
library(PopGenReport)
library(pegas)
library(diveRsity)

#################################################
#     Run genetic diversity analyses using sapply   #
#################################################
#create a species list
setwd("/home/user/git/Quercus_IUCN_samp_sims/")
species_list <- c("q_acerifolia",
                  "q_arkansana",
                  "q_austrina",
                  "q_boyntonii",
                  "q_carmenensis",
                  "q_cedrosensis",
                  "q_engelmannii",
                  "q_georgiana",
                  "q_graciliformis",
                  "q_havardii",
                  "q_hinckleyii",
                  "q_oglethorpensis",
                  "q_pacifica",
                  "q_tomentella")


# Function for calculating 3 required metrics (number of alleles, Ar, and He)
genAnalyze <- function(genindFile){
  # Read in genind file
  genind <- read.genepop(genindFile, ncode = 3)
  # Summarize the genind file and calculate the mean number of alleles
  sp_sum <- summary(genind)
  Nalleles <- mean(sp_sum$pop.n.all)
  # Calculate allelic richness
  Ar <- mean(colMeans(allelic.richness(genind)$Ar))
  #save mean for final output table
  # sp_gendiv_df[ind,3] <- mean(sp_poppr[1:length(levels(sp_genind@pop)),10])
  # Use Hs to calculate heterozygosity 
  # (Note that this generates slightly different values than the above call, by poppr)
  He <- mean(Hs(genind))
  return(c(Nalleles,Ar,He))
}

# Generate a list of the 2 folders of simulations results to parse
folder_list <- c("Alternative_sim_files", "Original_sim_files")

# Build a matrix that captures the mean results across replicates for each species (and for both folders)
allspecies_gendiv_df <- matrix(nrow = length(species_list), ncol = 6)

###loop over original and alternate simulations
for(folder in 1:length(folder_list)){
  ##loop over all 14 species
  for(sp in 1:length(species_list)){
    #move into each species folder
    setwd(paste0("/home/user/git/Quercus_IUCN_samp_sims/",folder_list[[folder]], "/", species_list[[sp]]))
    #load in genind objects
    genind_list = list.files(pattern = ".gen$")
    sp_gendiv_df <- matrix(nrow = length(genind_list), ncol = 3)
    
    # Use sapply (like lapply, but it "simplifies") to apply the genAnalyze function to the genind_list
    # For every filename in the list, 3 metrics are calculated. These are written to a matrix
    sp_gendiv_df <- t(sapply(genind_list, genAnalyze, USE.NAMES = FALSE))
    
    #save individual information
    if(folder == 1){
      
      allspecies_gendiv_df[sp,1:3] <- apply(sp_gendiv_df, 2, mean)
      
    }else{
      
      allspecies_gendiv_df[sp,4:6] <- apply(sp_gendiv_df, 2, mean)
    }
    write.csv(allspecies_gendiv_df, "allspecies_gendiv_df.csv")
  }
}
setwd("/home/user/git/Quercus_IUCN_samp_sims/")
write.csv(allspecies_gendiv_df, "allspecies_gendiv_df.csv")
