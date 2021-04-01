#Sampling script for all quercus species

#Library functions
library(adegenet)
library(car)
library(diveRsity)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(hierfstat)

#containing sub-folders
my_dir = "C:\\Users\\eschumacher\\Documents\\GitHub\\Quercus_IUCN_samp_sims\\Simulations"
setwd(my_dir)
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")

###function
source("G:/Shared drives/Emily_Schumacher/ten_oaks_gen/Fa_sample_funcs.R")
##functions
colMax <- function(data) sapply(data, max, na.rm = TRUE)
sample.pop<-function(genind_obj,vect_pop_ID,vect_samp_sizes){
  p<-length(vect_pop_ID)
  if (p>1) {
    for (p in 1:length(vect_pop_ID))
      alleles[p,]<-colSums(genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),],na.rm=T)
    alleles<-colSums(alleles)
  } else {alleles<-colSums(genind_obj[[vect_pop_ID[p]]]@tab[sample(1:nrow(genind_obj[[vect_pop_ID[p]]]@tab), vect_samp_sizes[p]),],na.rm=T)}
  
  alleles
}    
##list frequency category
list_allele_cat<-c("global","glob_v_com","glob_com","glob_lowfr","glob_rare","reg_rare","loc_com_d1","loc_com_d2","loc_rare")


#creating a list of the species we have simulated
species_list = c("\\q_acerifolia",
                 "\\q_arkansana",
                 "\\q_boyntonii",
                 "\\q_carmenensis",
                 "\\q_cedrosensis",
                 "\\q_engelmannii",
                 "\\q_georgiana",
                 "\\q_graciliformis",
                 "\\q_havardii",
                 "\\q_hinckleyii",
                 "\\q_oglethorpensis",
                 "\\q_pacifica")


##saving allelic richness
all_existing_by_sp_reps <- array(dim = c(100,length(list_allele_cat),length(species_list)))

##Data frames to store allele capture
all_cap_mean <- matrix(nrow = length(species_list), ncol = length(list_allele_cat))
all_cap_per <- matrix(nrow = length(species_list), ncol = length(list_allele_cat))

##mean results
all_mean_species <- matrix(nrow = length(species_list), ncol = length(list_allele_cat))
all_cap_mean <- matrix(nrow = length(species_list), ncol = length(list_allele_cat))
all_cap_per <- matrix(nrow = length(species_list), ncol = length(list_allele_cat))

##setting up output data frames
all_cap_df <- matrix(nrow = length(species_list), ncol = length(list_allele_cat))

##species name list 
for(a in 1:length(species_list)) {
  setwd(paste(mydir,species_list[a],sep=""))
  list_files = list.files(paste(mydir,species_list[a],sep=""), pattern = ".gen$")
  ##alleles existing by species for acerifolia 
  for(b in 1:100){
    
    #creating a temporary genind object (using Adegenet package) for each simulation replicate
    temp_genind = read.genepop(list_files[[b]], ncode=3) 
    
    n_ind <- table(temp_genind@pop)
    
    ##Then create a genpop file for temp_genind
    Spp_tot_genpop <- genind2genpop(temp_genind)
    
    ##separate by population 
    Spp_tot_genind_sep <- seppop(temp_genind)
    
    #
    all_cat <- get.allele.cat(Spp_tot_genpop, c(1:5), 2, n_ind)
    
    #defining the first and last individuals of the entire population, so we know where to sample between
    first_ind = 1
    last_ind = sum(table(temp_genind@pop)) 
    
    #for each replicate, sample up to 500 individuals, starting with 1
    for(c in 1:max_sample_size) {
      
      #this is a check to make sure that k, the sample size, doesn't exceed the species total pop. size
      #the loop will break if k is greater than the total pop. size
      #in other words, sampling stops once the entire population has been sampled, or when 500 samples is reached
      if(c <= sum(table(temp_genind@pop))) {
        #choosing which rows of the matrix to sample from
        #rows indicate individuals
        rows_to_samp = sample(first_ind:last_ind, c)
      }
    }
    #now calculate the alleles captured by sampling
    all_cap <- colSums(temp_genind@tab[rows_to_samp,], na.rm = T)
    
    for(d in 1:length(all_cat)) all_existing_by_sp_reps[b,d,a] <- sum((all_cat[[d]])>0,na.rm=T)
    
    ##all the high migration allele # per category
    for (e in 1:length(all_cat)) all_cat_cap[b,e,a] <- round(sum(all_cap[all_cat[[e]]]>0))
    
    for (f in 1:length(all_cat)) all_cat_per[b,f,a] <- round(sum(all_cap[all_cat[[f]]]>0)/length(all_cat[[f]]),4)
    
    
  }
}


##loops for taking mean results 
for(j in 1:length(all_cat)) {
  for(k in 1:12){
    
    all_mean_species[k,j] <- round(mean(all_existing_by_sp_reps[,j,k]), 3)
    all_cap_mean[k,j] <- round(mean(all_cat_cap[,j,k]), 3)
    all_cap_per[k,j] <- round(mean(all_cat_per[,j,k]),3)*100
    
  }
}

rownames(all_mean_species) <- species_list
colnames(all_mean_species) <- list_allele_cat

##loop to merge information
for(m in 1:length(species_list)){
  for(n in 1:length(all_cat)){
    all_cap_df[m,n] <- paste0(all_cap_per[m,n], "%", " ", "(", all_cap_mean[m,n], ")")
  }
}  

###add names and rows
rownames(all_cap_df) <- species_list
colnames(all_cap_df) <- list_allele_cat

write.csv(all_cap_df, "G:\\Shared drives\\Emily_Schumacher\\sampling_tenoaks\\allele_cap_12oaks.csv")
write.csv(all_mean_species, "G:\\Shared drives\\Emily_Schumacher\\sampling_tenoaks\\allele_cap_mean_12oaks.csv")
