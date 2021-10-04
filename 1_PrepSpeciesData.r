library(RSQLite)

# create an output directory if it doesn't exist
outputDirectory <- here::here("_data","species",sp_code)
ifelse(!dir.exists(outputDirectory), dir.create(outputDirectory), FALSE)

# get species information from the database
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
sp_data <- dbGetQuery(db_cem, paste("SELECT * FROM lkpSpecies WHERE CUTECODE = " , sQuote(sp_code), sep="") )
dbDisconnect(db_cem)



# get the training data ###############################################################################################
species <- read.table(here::here("_data","occurrence","Maxent_practiceinputspecies.csv"), header=TRUE,  sep=',')
species <- species[which(species$species==sp_code),]
species_pts <- species[,-1] # drop the 
