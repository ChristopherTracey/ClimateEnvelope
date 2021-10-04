


# get the training data ###############################################################################################
species <- read.table(here::here("_data","occurrence","Maxent_practiceinputspecies.csv"), header=TRUE,  sep=',')
species <- species[which(species$species==sp_code),]
species_pts <- species[,-1] # drop the 
