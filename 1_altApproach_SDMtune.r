
# create an output directory if it doesn't exist
ifelse(!dir.exists(here::here("_data","species",sp_code)), dir.create(here::here("_data","species",sp_code)), FALSE)

# get species information from the database
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
sp_data <- dbGetQuery(db_cem, paste("SELECT * FROM lkpSpecies WHERE CUTECODE = " , sQuote(sp_code), sep="") )
dbDisconnect(db_cem)

model_run_name <- paste0(sp_code, "_" , gsub(" ","_",gsub(c("-|:"),"",as.character(Sys.time()))))

# get the predictor data ##############################################################################################
cat("Loading the predictor data...")

# get the study area
studyArea <- arc.open(studyArea)
studyArea <- arc.select(studyArea)
studyArea <- arc.data2sf(studyArea)

# stack the predictors
predictors_Current <- stack(list.files(pathPredictorsCurrent, pattern = 'tif$', full.names=TRUE ))
predictors_Current <- projectRaster(predictors_Current, crs=crs(studyArea))
predictors_Future <- stack(list.files(pathPredictorsFuture, pattern = 'tif$', full.names=TRUE ))
predictors_Future <- projectRaster(predictors_Future, crs=crs(studyArea))

# check to see if the names are the same
setdiff(names(predictors_Current), names(predictors_Future))
setdiff(names(predictors_Future), names(predictors_Current))

# get set of random points for correlation analysis
set.seed(25)
bg <- dismo::randomPoints(predictors_Current, 10000)
bg <- prepareSWD(species="Bgs", a=bg, env=predictors_Current)
plotCor(bg, method="spearman", cor_th=0.7)
corVar(bg, method="spearman", cor_th=0.7)

# get the species data ################################################################################################
spData <- arc.open(spData_path)
spData <- arc.select(spData) #, dQuote(paste("SNAME=", sp_data$SNAME, sep=""))
spData <- arc.data2sf(spData)

spData_pro <- st_transform(spData, crs=crs(predictors_Current))

# write a shapeefile of the training data to the input folder for backup or other use
ifelse(!dir.exists(here::here("_data","species",sp_code,"input")), dir.create(here::here("_data","species",sp_code,"input")), FALSE)
ifelse(!dir.exists(here::here("_data","species",sp_code,"output")), dir.create(here::here("_data","species",sp_code,"output")), FALSE)
# write.csv(species_pts, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.csv")))
st_write(spData, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.shp")), append=FALSE)

# get some metadata


# convert species points to a lat/long df  MOVE THIS BELOW???
coords_pres <- data.frame(st_coordinates(spData_pro[,1]))
names(coords_pres) <- c("x","y")
md_ptTraining <- nrow(coords_pres) # get some metadata

# thin out the points to only one per cell
coords_pres_thin <- thinData(coords_pres, predictors_Current)
md_ptTrainingThinned <- nrow(spData) # get some metadata

# background coordinates
# These are 500 random locations, used as in place of absence values as 
# 'pseudoabsences' (the species probably doesn't occur at any random point)
coords_bg <- as.data.frame(dismo::randomPoints(predictors_Current, 500))
md_bg <- nrow(coords_bg) # get some metadata
coords_bg_thin <- thinData(coords_bg, predictors_Current)
md_bgThinned <- nrow(coords_bg_thin) # get some metadata


# # check to see if the coordinate systems are the same
# a <- stringr::str_split(as.character(crs(predictors_Current)), ' ')
# b <- stringr::str_split(as.character(crs(spData)), ' ')
# identical(sort(unlist(a)), sort(unlist(b)))
# rm(a,b)

# create SDW object ############################
data.SWD <- prepareSWD(species=unique(spData$SNAME), p=coords_pres, a=coords_bg, env=predictors_Current)
data.SWD

plotCor(data.SWD , method="spearman", cor_th=0.7)

################################################3333333333333333333
# run the models ##################################
#i = 2  # temp just for testing

for(i in 2:length(ModelMethods)){
  cat("--------------------------------------------------------\n")
  cat(paste("Running a ", ModelMethods[i], " model for ", unique(spData$SNAME),".\n", sep=""))
  # create folds for crossvalidated model
  folds <- randomFolds(data.SWD, k=4, only_presence=TRUE) #, seed=5
  cv_model <- train(ModelMethods[i], data=data.SWD)
  cv_model <- train(ModelMethods[i], data=data.SWD, folds=folds)
  cv_model
  auc(cv_model)
  auc(cv_model, test=TRUE)
  
  #model <- train(method=ModelMethods[i], data=data.SWD)  # model <- train(method = "Maxent", data = data, fc = "lh", reg = 0.5, iter = 700)
  #model
  
  # calculating 
  cat("- variable importance\n")
  if(ModelMethods[i]=="Maxent"){
  # vi <- maxentVarImp(model)
    # vi
    # plotVarImp(vi[, 1:2])
    #reduced_variables_model <- varSel(cv_model, metric="tss", bg4cor=bg, method="spearman", cor_th=0.7, use_pc=TRUE)
  } else if(ModelMethods[i]=="RF"|ModelMethods[i]=="BRT") {
    vi <- SDMtune::varImp(cv_model, permut = 10) #,
    vi
    plotVarImp(vi[, 1:2])
    
    # jacknifing approach to variable selection
    md_tssPre <- tss(cv_model, test=TRUE)
    md_aucPre <- auc(cv_model, test=TRUE)
    cat("Testing TSS before: ", md_tssPre, "\n")
    
   # reduced_variables_model <- varSel(cv_model, metric="tss", bg4cor=bg, method="spearman", cor_th=0.7)
    reduced_variables_model <- reduceVar(cv_model, th=5, metric="tss", permut=1) #, use_jk=TRUE
    reduced_variables_model <- reduceVar(cv_model, metric="auc", th=10, use_jk = TRUE)
    
    
    cat("The following variables were used in the final model:", names(reduced_variables_model@data@data), "\n")
    md_tssPost <- tss(reduced_variables_model, test=TRUE)
    md_aucPost <- auc(reduced_variables_model, test=TRUE)
    cat("Testing TSS after: ", md_tssPost, "\n")
    
  } else {
    cat("No valid variable importance method exists...")
  }
  
  # insert some model run metadata
  sf_metadata <- data.frame("sp_code"=sp_code, "model_run_name"=model_run_name, "model_type"=ModelMethods[i], "modeller"=modeller, "TrainingPoints"=md_ptTraining, "TrainingPoints_thinned"=md_ptTrainingThinned, "BackgroundPoints"=md_bg, "BackgroundPoints_thinned"=md_bgThinned, "AUCpre"=md_aucPre, "AUCpost"=md_aucPost, "TSSpre"=md_tssPre, "TSSpost"=md_tssPost)  
  
  db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
  SQLquery <- paste("INSERT INTO model_runs (", paste(names(sf_metadata), collapse = ',') ,") VALUES (",paste(sQuote(sf_metadata[1,]), collapse = ','),");", sep="") #sp_code, model_run_name, modeller,TrainingPoints
  dbExecute(db_cem, SQLquery )
  dbDisconnect(db_cem)
  
  # predict the model to the current env predictors
  cat("- predicting the model to the current env predictors\n")
  timeframe <- "current"
  if(ModelMethods[i]=="Maxent"){
    map <- predict(reduced_variables_model, data=predictors_Current, type="cloglog")
  } else if(ModelMethods[i]=="RF"|ModelMethods[i]=="BRT") {
    map <- predict(reduced_variables_model, data=predictors_Current)
  } else {
    cat("No valid model predictor method exists...")
  }
  plotPred(map) # plot and write the current map
  writeRaster(map, here::here("_data","species",sp_code,"output",paste(model_run_name, "_", ModelMethods[i], "_", timeframe, ".tif", sep="")),  "GTiff", overwrite=TRUE)
  
  # predict the model to the future env predictors
  cat("- predicting the model to the future env predictors\n")
  timeframe <- "future"
  if(ModelMethods[i]=="Maxent"){
    map <- predict(reduced_variables_model, data=predictors_Future, type="cloglog")
  } else if(ModelMethods[i]=="RF"|ModelMethods[i]=="BRT") {
    map <- predict(reduced_variables_model, data=predictors_Future)
  } else {
    cat("No valid model predictor method exists...")
  }
  plotPred(map) # plot and write the future map
  writeRaster(map, here::here("_data","species",sp_code,"output",paste(model_run_name, "_", ModelMethods[i], "_", timeframe, ".tif", sep="")),  "GTiff", overwrite=TRUE)
  
  # cleanup
  #rm(cv_model, reduced_variables_model)
  cat("- finished with the model\n")
}
