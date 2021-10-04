library(RSQLite)

# create an output directory if it doesn't exist
ifelse(!dir.exists(here::here("_data","species",sp_code)), dir.create(here::here("_data","species",sp_code)), FALSE)

# get species information from the database
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
sp_data <- dbGetQuery(db_cem, paste("SELECT * FROM lkpSpecies WHERE CUTECODE = " , sQuote(sp_code), sep="") )
dbDisconnect(db_cem)


# get the training data ###############################################################################################
species <- read.table(here::here("_data","occurrence","Maxent_practiceinputspecies.csv"), header=TRUE,  sep=',')
species <- species[which(species$species==sp_data$SNAME),]
species_pts <- species[,-1] # drop the 
rm(species)

# write a csv of the training data to the input folder for backup or other use
ifelse(!dir.exists(here::here("_data","species",sp_code,"input")), dir.create(here::here("_data","species",sp_code,"input")), FALSE)
write.csv(species_pts, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.csv")))

