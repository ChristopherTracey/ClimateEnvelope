


# get the training data ###############################################################################################
species <- read.table(here::here("Maxent_practiceinputspecies.csv"), header=TRUE,  sep=',')
species <- species[which(species$species==sp_code),]
species_pts <- species[,-1] # drop the 
