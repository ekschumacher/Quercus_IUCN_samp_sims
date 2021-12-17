##create a list to store all occurrence records 
oak_occ <- list()

##list out all species 
oak_names <- c("QUAC", "QUAR", "QUAU", "QUBO",
               "QUCA", "QUCE", "QUEN", "QUGE",
               "QUGR", "QUHA", "QUHI", "QUOG",
               "QUPA", "QUTO")

##project spdf list 
oak_pcs_spdf_list <- list()

##list for saving occurrences 
occurence_no_auto <- list()

##list of data frames 
oak_red_dfs <- list()

##Projection string
projection <- c("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

for(sp in 1:length(oak_names)){
  
  ##load in data frames for presence 
  oak_occ <- read.csv(paste0(oak_names[[sp]], "/", oak_names[[sp]], "_occ.csv"))
  
  ##create point data frame
  oak_occ_spdf <- matrix(nrow = length(oak_occ[,1]), ncol = length(oak_occ[,1]))
  
  for(lon in 1:length(oak_occ[,1])){
    for(lat in 1:length(oak_occ[,1])){
      
      oak_occ_spdf[lon,lat] <- (distm(oak_occ[lon,2:3], oak_occ[lat,2:3])/1000) < 50
      
    }
  }
  
  ##create tri matrix with NAs for bottom rows
  oak_occ_spdf[lower.tri(oak_occ_spdf, diag=TRUE)] <- NA
  
  ##now determine where points within 1 km are 
  auto_cols <- colSums(oak_occ_spdf, na.rm=TRUE) == 0
  
  ##now subset the original data by this 
  occurence_no_auto <- oak_occ[auto_cols,]
  
  ##create data frames of reduced individual data frames 
  oak_red_dfs <- data.frame(occurence_no_auto)
  
  ##write out 
  write.csv(oak_red_dfs, paste0(oak_names[[sp]], "/", oak_names[[sp]], "_occ_red.csv"))
  
}
