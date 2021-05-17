# all_quercus_sampling.R - Sampling script for all quercus species

#First, this script imports simulation files in .arp format and converts them to .gen format.
#Then, the script runs the main sampling loop. This loop runs through every replicate of each species and
#samples from 1 to 500 individuals for each simulation replicate. From this, we calculate the proportion of
#alleles captured by the sample size, by dividing the captured alleles by the total alleles present in the
#simulation. The results are saved in an array, and then saved as .Rdata, which can be imported to the next
#R script, all_quercus_processing
#The loop also implements Fst calculations, which can be turned off using a flag: fst_flag
#similarly, the file conversion function can be turned off once files have been converted to .gen once
#(there is no need to convert more than once, unless simulation files are changed) using the flag: conversion_flag

#This script was written in collaboration by Kaylee Rosenberger, Emily Schumacher, and Dr. Sean Hoban

########################################################################################################
#Library functions
library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(hierfstat)

#FLAGS
#file conversion flag
#set to true once files have been converted once
#false if you want to convert files
conversion_flag = FALSE
#Fst flag
#Fst code adds a lot of time to run the code 
#so if you don't want to run it, keep Fst off by setting it FALSE
fst_flag = TRUE

#Set working directory
mydir = "C:\\Users\\kayle\\Documents\\Quercus_IUCN_samp_sims_local\\Simulation_files"
setwd(mydir)

#creating a list of the species we have simulated
species_list = c("\\q_acerifolia",
                 "\\q_arkansana",
                 "\\q_austrina",
                 "\\q_boyntonii",
                 "\\q_carmenensis",
                 "\\q_cedrosensis",
                 "\\q_engelmannii",
                 "\\q_georgiana",
                 "\\q_graciliformis",
                 "\\q_havardii",
                 "\\q_hinckleyii",
                 "\\q_oglethorpensis",
                 "\\q_pacifica",
                 "\\q_tomentella")

##allele frequency categories 
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

#defining the maximum number of individuals we want to sample
#for practical purposes, this will be 500 indivduals 
max_sample_size = 500

#number of replicates of genetic simulation
num_replicates = 1000

#Defining an import function
#converts all arlequin files in a folder to genepop files
import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath)
  temp_list_1 = list.files(mypath, mypattern)
  temp_list_2 = list(length = length(temp_list_1))
  for(z in 1:length(temp_list_1)){temp_list_2[[z]]=arp2gen(temp_list_1[z])}
  temp_list_2
}

#converting all simulation files from arlequin format to genepop format using defined import function
if(conversion_flag == TRUE) {
  for(x in 1:length(species_list)) {
    import_arp2gen_files(paste(mydir,species_list[x],sep=""),".arp$")
  }
}

#pre-defining the array to store results
#first dimension: 500, sampling from 1 to 500 individuals per species, saving results for each iteration
#second dimension: 100 for 100 simulation replicates per species
#third dimension: 14, for 14 quercus species. this is represented by the outer for loop (14 matrix 'slices')
final_quercus_results = array(0, dim = c(500,1000,14))

final_alleles_all_quercus = array(0, dim=c(500,1000,14))

#storing Fst results
#saving a list of genind objects created
temp_genind_list <- list()

#list of hierfstat
temp_hierfstat <- list()

##min, max, mean of replicates 
mean_max_min_fst = array(dim = c(3,100,14))#Fst run on 100 replicates

##
all_cap_samp <- array(dim = c(100, length(list_allele_cat), length(species_list)))

##
all_cap_samp_per <- array(dim = c(100, length(list_allele_cat), length(species_list)))

##species abbrevs 
species_abbrevs <- c("QUAC","QUAR","QUBO","QUCA","QUCE","QUEN","QUGE","QUGR","QUHA","QUHI","QUOG","QUPA")

##allelic existing by species table 
all_existing_by_sp_df <- matrix(nrow = length(species_list), ncol = length(list_allele_cat))
colnames(all_existing_by_sp_df) <- list_allele_cat
rownames(all_existing_by_sp_df) <- species_abbrevs

##alleles sampled 
all_cap_cat_per <- matrix(nrow = length(species_list), ncol = length(list_allele_cat))
colnames(all_cap_cat_per) <- list_allele_cat
rownames(all_cap_cat_per) <- species_abbrevs

##combo data frame with alleles captured by category vs. % captured 
all_cap_cat_df <- matrix(nrow = length(species_list), ncol = length(list_allele_cat))

###############################################################################################
#SAMPLING/Fst

#Loop to simulate sampling
#First, create a list of all genepop files (all replicates) to loop over
#the variable 'i' represents each replicate
for(i in 1:length(species_list)) {
  setwd(paste(mydir,species_list[i],sep=""))
  list_files = list.files(paste(mydir,species_list[i],sep=""), pattern = ".gen$")
  for(j in 1:100) { #(length(list_files))
    #creating a temporary genind object (using Adegenet package) for each simulation replicate
    temp_genind = read.genepop(list_files[[j]], ncode=3) 
   
       #First, calculate number of individuals per population
       n_ind <- table(temp_genind@pop)
       
       ##Then create a genpop file for temp_genind
       Spp_tot_genpop <- genind2genpop(temp_genind)
       
       # ##separate by population 
       Spp_tot_genind_sep <- seppop(temp_genind)
       
       ##get categories for all alleles captured 
       allele_cat_tot <- get.allele.cat(Spp_tot_genpop, c(1:5), 2, n_ind)
       
       ##calculate the total number of alleles in each frequency category over 9 allele categories
       for (a in 1:length(allele_cat_tot)) all_existing_by_sp_reps[j,a,i] <- sum(allele_cat_tot[[a]]>0,na.rm=T)
       
       ##create output for all alleles captured by 
       for(b in 1:length(allele_cat_tot)) all_existing_by_sp_df[i,b] <- mean(all_existing_by_sp_reps[,b,i])
       
      ###now start sampling code 
      #defining the first and last individuals of the entire population, so we know where to sample between
      first_ind = 1
      last_ind = sum(table(temp_genind@pop))
      
      #for each replicate, sample up to 500 individuals, starting with 1
       for(k in 1:max_sample_size) {
      # 
      #   #this is a check to make sure that k, the sample size, doesn't exceed the species total pop. size
      #   #the loop will break if k is greater than the total pop. size
      #   #in other words, sampling stops once the entire population has been sampled, or when 500 samples is reached
         if(k <= sum(table(temp_genind@pop))) {
           #choosing which rows of the matrix to sample from
           #rows indicate individuals
           rows_to_samp = sample(first_ind:last_ind, k)
       
           #saving the alleles sampled
           if(k == 1) {
             sample_n_alleles = sum(temp_genind@tab[rows_to_samp,]>0)
           } else {
             sample_n_alleles = sum(colSums(temp_genind@tab[rows_to_samp,])>0)
           }
       
      #     #calculating the total alleles
      #     total_alleles = ncol(temp_genind@tab)
      # 
      #     #saving the proportion of alleles captured -> alleles sampled/total alleles
      #     #represents genetic conservation success
      #     final_quercus_results[k,j,i] = sample_n_alleles/total_alleles 
      # 
      #     #saving the total alleles present across the populations for each species, and each replicate
      #     final_alleles_all_quercus[k,j,i] = total_alleles 
           
       }else {
        break
       }
         
       }
      ##start adding allelic capture code 
      alleles_cap <- colSums(temp_genind@tab[rows_to_samp,], na.rm = T)
      
      ##allelic capture code 
      for (b in 1:length(allele_cat_tot)) all_cap_samp[j,b,i] <- sum(alleles_cap[allele_cat_tot[[b]]]>0)
      
      ##percent captured data frame 
      all_cap_samp_per[,,i] <- all_cap_samp[,,i]/all_existing_by_sp_reps[,,i]
      
      ##all_cap_cat
      for (c in 1:length(allele_cat_tot)) all_cap_cat_per[i,c] <- round(mean(all_cap_samp_per[,c,i])*100,3)
      
      ##loop to write out nice df 
      for(d in 1:length(species_abbrevs)){
        for(e in 1:length(allele_cat_tot)){
          all_cap_cat_df[d,e] <- paste0(signif(all_cap_cat_per[d,e], 3), "%", " ", "(", signif(all_existing_by_sp_df[d,e],3), ")")
        }
      }  
      
       
   }
}

##write out alleles existing within each categories  
setwd("C:\\Users\\eschumacher\\Documents\\GitHub\\Quercus_IUCN_samp_sims\\Results")
write.csv()

