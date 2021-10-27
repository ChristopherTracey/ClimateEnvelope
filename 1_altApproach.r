

# create an output directory if it doesn't exist
ifelse(!dir.exists(here::here("_data","species",sp_code)), dir.create(here::here("_data","species",sp_code)), FALSE)

# get species information from the database
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
sp_data <- dbGetQuery(db_cem, paste("SELECT * FROM lkpSpecies WHERE CUTECODE = " , sQuote(sp_code), sep="") )
dbDisconnect(db_cem)

model_run_name <- paste0(sp_code, "_" , gsub(" ","_",gsub(c("-|:"),"",as.character(Sys.time()))))



# get the species data ################################################################################################
spData <- arc.open(spData)
spData <- arc.select(spData)
spData <- arc.data2sf(spData)

# write a shapeefile of the training data to the input folder for backup or other use
ifelse(!dir.exists(here::here("_data","species",sp_code,"input")), dir.create(here::here("_data","species",sp_code,"input")), FALSE)
ifelse(!dir.exists(here::here("_data","species",sp_code,"output")), dir.create(here::here("_data","species",sp_code,"output")), FALSE)
# write.csv(species_pts, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.csv")))
st_write(spData, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.shp")), append=FALSE)

# get some metadata
md_ptTraining <- nrow(spData)


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
sp_coords <- data.frame(st_coordinates(spData[,1]))
names(sp_coords) <- c("lon","lat")


# get the predictor data ##############################################################################################
cat("Loading the predictor data...")
predictors_Current <- stack(list.files(pathPredictorsCurrent, pattern = 'tif$', full.names=TRUE ))
predictors_Current <- projectRaster(predictors_Current, crs=crs(studyArea))
predictors_Future <- stack(list.files(pathPredictorsFuture, pattern = 'tif$', full.names=TRUE ))
predictors_Future <- projectRaster(predictors_Future, crs=crs(studyArea))

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




# PREP THE DATA for the SDM package ##############################################

# get a list of variables
VarList <- names(predictors_Current)
VarListModelEq <- paste(VarList, collapse=" + ")

# Prep data format for the sdm package
sdm.pkg.df_pres <- cbind(sp_coords, presVals)
sdm.pkg.df_pres$Ey <- 1
names(sdm.pkg.df_pres)[1:2] <- c("x", "y")
sdm.pkg.df_abs <- data.frame(cbind(backgr, absVals))
sdm.pkg.df_abs$Ey <- 0
sdmdf_sdmpkg <- rbind(sdm.pkg.df_pres, sdm.pkg.df_abs)
eq <- as.formula(paste0("Ey ~ ", VarListModelEq))
sdmdata_sdmpkg <- sdmData(eq, train=sdmdf_sdmpkg)

# get shared data




#also get correlated env var information
db_cem <- dbConnect(SQLite(), dbname=nm_db_file)
SQLquery <- "SELECT rasName, CorrGroup FROM lkpEnvVar WHERE CorrGroup IS NOT NULL order by CorrGroup;"
corrdEVs <- dbGetQuery(db_cem, statement = SQLquery)
dbDisconnect(db_cem)
rm(db_cem, SQLquery)


b <- sdm(eq, data=sdmdata_sdmpkg, methods=c('rf','glm','maxent'), replicatin='sub', test.percent=30, n=2)

b
roc(b, smooth = TRUE)
b.var <- getVarImp(b)

#####################################################################################
# General Linear Model


# Run a GLM model using the sdm package
sdm_ml.glm <- sdm::sdm(eq, data=sdmdata_sdmpkg, methods=c("glm"))
prediction_ml.glm <- raster::predict(sdm_ml.glm, predictors_Current)
project.sdm(prediction_ml.glm, "GLM SDM")

rcurve(sdm_ml.glm)
roc(sdm_ml.glm)
a <- getVarImp(sdm_ml.glm)

a@varImportance


# predict to the future
prediction_ml.glm_future <- raster::predict(sdm_ml.glm, predictors_Future)
project.sdm(prediction_ml.glm_future, "GLM SDM Future")

#####################################################################################
# RandomForest
# Now Let's try more strictly Machine Learning: Random Forest
sdm::installAll()
library(sdm)

# Run the model and project
sdm_rf <- sdm::sdm(eq, data=sdmdata_sdmpkg, methods=c("rf"))
sdm_rf <- sdm::sdm(eq, data=sdmdata_sdmpkg, methods=c("rf"), replication='sub', test.percent=30, n=10)
prediction_rf <- raster::predict(sdm_rf, predictors_Current)
project.sdm(prediction_rf, "Random Forest SDM")
rf_VarImport <- getVarImp(sdm_rf)


prediction_rf_future <- raster::predict(sdm_rf, predictors_Future)
project.sdm(prediction_rf_future, "Random Forest Future SDM")

#####################################################################################
# Maximum Entropy

# Maxent, need to install maxent (https://biodiversityinformatics.amnh.org/open_source/maxent/)
# and place it here:
system.file("java", package="dismo")

cat("Running a test model for variable filtering...")
cem_MaxEntTest <- maxent(predictors_Current, sp_coords, silent=TRUE)
# remove the least important variables
me.dat <- as.data.frame(slot(cem_MaxEntTest, "results"))
me.imp.dat <- me.dat[grepl("permutation.importance",rownames(me.dat)), ,drop = FALSE]
me.imp.dat <- cbind(me.imp.dat, "var" = unlist(lapply(rownames(me.imp.dat), FUN = function(x) strsplit(x, "\\.")[[1]][[1]])))

impvals <- me.imp.dat
names(impvals) <- c("imp","var")

OriginalNumberOfEnvars <- nrow(impvals)

# correlated variables are gotten above for all models
if(nrow(corrdEVs) > 0 ){
  for(grp in unique(corrdEVs$CorrGroup)){
    vars <- corrdEVs[corrdEVs$CorrGroup == grp,"rasName"]
    imp.sub <- impvals[impvals$var %in% vars,, drop = FALSE]
    varsToDrop <- imp.sub[!imp.sub$imp == max(imp.sub$imp),, drop = FALSE]#suppressWarnings()
    impvals <- impvals[!impvals$var %in% varsToDrop$var,,drop = FALSE]
  }
  rm(vars, imp.sub, varsToDrop)
}

# set the percentile, here choosing above 25% percentile
envarPctile <- 0.25
y <- quantile(impvals$imp, probs = envarPctile)
impEnvVars <- impvals[impvals$imp > y,]
subsetNumberofEnvars <- nrow(impEnvVars)
rm(y)
# which columns are these, then flip the non-envars to TRUE
impEnvVarCols <-names(predictors_Current)[which(names(predictors_Current) %in% impEnvVars$var)]
# subset!
predictors_Current1 <- raster::subset(predictors_Current, impEnvVarCols)
predictors_Future1 <- raster::subset(predictors_Future, impEnvVarCols)

# reset the indvarcols object
#indVarCols <- c(6:length(names(me.df.full)))
rm(impvals, impEnvVars, impEnvVarCols)

# ########################
sdm_maxent <- maxent(predictors_Current1, sp_coords)
prediction_maxent <- dismo::predict(sdm_maxent, predictors_Current1)
project.sdm(prediction_maxent, "MaxEnt SDM ")

prediction_maxent_future <- dismo::predict(sdm_maxent, predictors_Future1)
project.sdm(prediction_maxent_future, "MaxEnt SDM Future")

# Look at response for each predictor
response(sdm_maxent)

# Look at variable contribution for maxent
plot(sdm_maxent)


#############################################
# insert metadata into database
sf_metadata <- data.frame("sp_code"=sp_code, "model_run_name"=model_run_name, "modeller"=modeller, "md_ptTraining"=md_ptTraining)  
  db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
  SQLquery <- paste("INSERT INTO model_runs (sp_code, model_run_name, modeller,TrainingPoints) VALUES (",paste(sQuote(sf_metadata[1,]), collapse = ','),");", sep="")
  dbExecute(db_cem, SQLquery )
  dbDisconnect(db_cem)

rm(md_ptTraining)

