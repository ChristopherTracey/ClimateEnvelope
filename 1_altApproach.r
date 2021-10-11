

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


# get the species data ################################################################################################
spData <- arc.open(spData)
spData <- arc.select(spData)
spData <- arc.data2sf(spData)

# write a csv of the training data to the input folder for backup or other use
# ifelse(!dir.exists(here::here("_data","species",sp_code,"input")), dir.create(here::here("_data","species",sp_code,"input")), FALSE)
# ifelse(!dir.exists(here::here("_data","species",sp_code,"output")), dir.create(here::here("_data","species",sp_code,"output")), FALSE)
# write.csv(species_pts, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.csv")))

# create a sf points layer
# species_sf <- st_as_sf(species_pts, coords=c("lon","lat")) # 
# plot(species_sf)
# species_sf <- st_transform(species_sf, projPA)

st_write(spData, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.shp")), append=FALSE)


# get the basemap data ################################################################################################
studyArea <- arc.open(studyArea)
studyArea <- arc.select(studyArea)
studyArea <- arc.data2sf(studyArea)

# check to see if the coordinate systems are the same
a <- stringr::str_split(as.character(crs(studyArea)), ' ')
b <- stringr::str_split(as.character(crs(spData)), ' ')
identical(sort(unlist(a)), sort(unlist(b)))
rm(a,b)
# Let's take a look at our raw occurrence points
ggplot() +
  geom_sf(studyArea, mapping=aes(fill=Id)) +
  geom_sf(spData, mapping=aes(col=EORANK))


# convert species points to a lat/long df  MOVE THIS BELOW???
sp_coords < data.frame(st_coordinates(spData[,1]))
names(sp_coords) <- c("lon","lat")


# get the predictor data ##############################################################################################
cat("Loading the predictor data...")
predictors_Current <- stack(list.files(here::here("_data","env_vars","Climate_Current"), pattern = 'asc$', full.names=TRUE ))
predictors_Future <- stack(list.files(here::here("_data","env_vars","Climate_Future2050rfp45"), pattern = 'asc$', full.names=TRUE ))

# check to see if the names are the same
setdiff(names(predictors_Current), names(predictors_Future))
setdiff(names(predictors_Future), names(predictors_Current))

# extract values to point from the raster stack
presVals <- raster::extract(predictors_Current, spData, method="simple")

# These are 500 random locations, used as in place of absence values as 
# 'pseudoabsences' (the species probably doesn't occur at any random point)
backgr <- randomPoints(predictors_Current, 500)
# predictor values at random locations
absVals <- raster::extract(predictors_Current, backgr)

# We know that expected habitat suitability (Ey) is 1 for areas where the species
# was found, and we assume it's 0 for the random background points
Ey <- c(rep(1, nrow(presVals)), rep(0, nrow(absVals)))

# Now here we have a dataframe with the response variable (Ey) and corresponding
# predictor values
sdmdata <- data.frame(cbind(Ey, rbind(presVals, absVals)))
View(sdmdata)


# Let's start Modeling! Here's the most basic SDM

## Linear Model SDM
# E(y) = mx + b + e
# E(y) is expected habitat suitability
# b = intercept
# m = coefficient
# x = predictor
# e = error

# Assumes: Normal distributed error
#          mean = 0
#          Variance(y) is constant



# Let's make  GLM with all 3 variables
par(mfrow=c(1,1))
sdm_glm3 <- glm(Ey ~ CMD + DD5 + MAP + MAR + MSP + PAS + RH + TD, data=sdmdata, family = binomial)
prediction_glm3 <- raster::predict(predictors_Current, sdm_glm3)
prediction_glm3.Ey <- exp(prediction_glm3) / (1+exp(prediction_glm3))
project.sdm(prediction_glm3.Ey, "GLM SDM (D. californica)")


# RandomForest
# Now Let's try more strictly Machine Learning: Random Forest
sdm::installAll()
library(sdm)

# Prep data format for the sdm package
sdm.pkg.df_pres <- cbind(sp_coords, presVals)
sdm.pkg.df_pres$Ey <- 1
names(sdm.pkg.df_pres)[1:2] <- c("x", "y")
sdm.pkg.df_abs <- data.frame(cbind(backgr, absVals))
sdm.pkg.df_abs$Ey <- 0
sdmdf_sdmpkg <- rbind(sdm.pkg.df_pres, sdm.pkg.df_abs)
sdmdata_sdmpkg <- sdmData(Ey ~ CMD + DD5 + MAP + MAR + MSP + PAS + RH + TD, train = sdmdf_sdmpkg)


# Run the model and project
sdm_rf <- sdm::sdm(Ey ~ CMD + DD5 + MAP + MAR + MSP + PAS + RH + TD, data=sdmdata_sdmpkg, methods=c("rf"))
prediction_rf <- raster::predict(sdm_rf, predictors_Current)
project.sdm(prediction_rf, "Random Forest SDM (D. californica)")
getVarImp(sdm_rf)

prediction_rf_future <- raster::predict(sdm_rf, predictors_Future)
project.sdm(prediction_rf_future, "Random Forest Futue SDM (D. californica)")

