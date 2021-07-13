##########################
######## Libraries #######
##########################
library(adegenet)
library(hierfstat)

#####################################
############ Directories ############
#####################################
##data file input paths
allpop_genind_path <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Ten_Oaks\\DataFiles\\genind_files\\allpop\\"
tenoak_df_path <- "G:\\My Drive\\Hoban_Lab_Docs\\Projects\\Ten_Oaks\\DataFiles\\Oak_Score_Df"

##species list 
oak_species_list <- c("QUAC","QUAJ","QUHI","QUPA")

#############################################
############ Conversion Code ################
#############################################
##all pops conversion 
setwd("allpop_genind_path")

##write in genind files 
allpop_arp_list <- list.files(pattern = "_wild_pop.arp$")

##convert

for(o in 1:length(allpop_arp_list)){
  
  arp2gen(allpop_arp_list[[o]])
  
}

#######################################
############ load in files ############
#######################################
setwd("C:\\Users\\eschumacher\\Documents\\GitHub\\Quercus_IUCN_samp_sims\\oak_geninds")

##write in genind files 
ten_oaks_genind_wild_list <- list.files(pattern = "_wild_pop.gen$")

##conversion list
tenoak_wild_genind_list <- list()

##read in data frame 
#tenoak_wild_df_list <- list.files(path = tenoak_df_path, 
    #                              pattern = "_wild_pops.csv$")

##storage list for data frames 
#tenoak_wild_dfs <- list()

##create data frame
fst_df <- matrix(nrow = 3, ncol = 7)

##loop to load in genind files
for(o in 1:length(ten_oaks_genind_wild_list)){
  
  ##load in genind files
  tenoak_wild_genind_list[[o]] <- read.genepop(ten_oaks_genind_wild_list[[o]], ncode = 3)
 
  #First put into poppr format
  popr_test <- as.genclone(tenoak_wild_genind_list[[o]])
  strata(popr_test) <- other(popr_test)$population_hierarchy[-1]
  list_a<-mlg.id(popr_test)
  #Function to pull out individual indices where clone length greater than 1
  clone_index<-which(sapply(list_a,function(x) length(x)>1))
  
  #This removes clones and then saves as new file for Genealex if desired
  popr_nocl<-clonecorrect(popr_test,strata=~Pop)
  #genind2genalex(genclone2genind(popr_nocl),file="QH_clone_free.csv")
  
  #Create genpop and genind objects that now have no clones- GI_nocl, GP_nocl
  GI_nocl<-genclone2genind(popr_nocl); 	GP_nocl<-genind2genpop(GI_nocl)
  
  fst_conversion <- genind2hierfstat(GI_nocl)
  
  ##calculate fst 
  fst_spp <- pairwise.neifst(fst_conversion)

  ##store in df
  fst_df[1,o] <- mean(sapply(fst_spp, mean, na.rm = T), na.rm = T)
  fst_df[2,o] <- min(fst_spp, na.rm = T)
  fst_df[3,o] <- max(fst_spp, na.rm = T)
}
