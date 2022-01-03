
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

# load and stack the Current predictors
predictors_Current <- stack(list.files(pathPredictorsCurrent, pattern = 'tif$', full.names=TRUE ))

# check to see if the crs of the current predictors and the study area are the same
a <- str_split(crs(studyArea), ' ')
b <- str_split(crs(predictors_Current), ' ')
if(isFALSE(identical(sort(unlist(a)), sort(unlist(b))))){
  cat("CRS's of the study area and predictors are mismatched. Reprojecting...")
  studyArea <- st_transform(studyArea, crs=crs(predictors_Current)) # reproject the study area
  a <- str_split(crs(studyArea), ' ')
  identical(sort(unlist(a)), sort(unlist(b)))
} else {
  cat("CRS's of the study area and predictors are the same. No worries...")
}
rm(a,b)

# load and stack the Future predictors
predictors_Future4.5 <- stack(list.files(pathPredictorsFuture4.5, pattern = 'tif$', full.names=TRUE ))
crs(predictors_Future4.5) <- crs(predictors_Current)
predictors_Future8.5 <- stack(list.files(pathPredictorsFuture8.5, pattern = 'tif$', full.names=TRUE ))
crs(predictors_Future8.5) <- crs(predictors_Current)

# check to see if the names are the same
if(isTRUE(identical(sort(names(predictors_Current)), sort(names(predictors_Future4.5))))&isTRUE(identical(sort(names(predictors_Current)), sort(names(predictors_Future8.5))))){
  cat("All raster layers are correctly named...")
} else {
  cat("You have a name mismatch somewhere in the future raster layers")
}

# get set of random points for correlation analysis
set.seed(25)
bg <- dismo::randomPoints(predictors_Current, 10000)
bg <- prepareSWD(species="Bgs", a=bg, env=predictors_Current)
plotCor(bg, method="spearman", cor_th=0.7)
corVar(bg, method="spearman", cor_th=0.7)

# get the species data ################################################################################################
spData <- arc.open(spData_path)
spData <- arc.select(spData)
spData <- arc.data2sf(spData)
SpData <- spData[spData$SNAME==sp_data$SNAME,]
SpData <- SpData[SpData$USEDATA %in% c("Y","y"),]#keep only the species marked as a "Y" for modeling

# check to see if projections are the same and modify if not
a <- str_split(crs(SpData), ' ')
b <- str_split(crs(predictors_Current), ' ')
if(isFALSE(identical(sort(unlist(a)), sort(unlist(b))))){
  cat("CRS's of the species data and predictors are mismatched. Reprojecting the species data...")
  spData_pro <- st_transform(SpData, crs=crs(predictors_Current)) # reproject the study area
  a <- str_split(crs(spData_pro), ' ')
  identical(sort(unlist(a)), sort(unlist(b)))
} else {
  cat("CRS's of the species data and predictors are the same. No worries...")
}
rm(a,b,spData,SpData)



# write a shapefile of the training data to the input folder for backup or other use
ifelse(!dir.exists(here::here("_data","species",sp_code,"input")), dir.create(here::here("_data","species",sp_code,"input")), FALSE)
ifelse(!dir.exists(here::here("_data","species",sp_code,"output")), dir.create(here::here("_data","species",sp_code,"output")), FALSE)
st_write(spData_pro, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.shp")), append=FALSE)

# convert species points to a lat/long df  MOVE THIS BELOW???
coords_pres <- data.frame(st_coordinates(spData_pro[,1]))
names(coords_pres) <- c("x","y")
md_ptTraining <- nrow(coords_pres) # get some metadata

# thin out the points to only one per cell
coords_pres_thin <- thinData(coords_pres, predictors_Current)
md_ptTrainingThinned <- nrow(spData_pro) # get some metadata

# background coordinates. These are 500 random locations, used as in place of absence values as 'pseudoabsences' (the species probably doesn't occur at any random point)
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
data.SWD <- prepareSWD(species=unique(spData_pro$SNAME), p=coords_pres_thin, a=coords_bg_thin, env=predictors_Current)
data.SWD

#in case things aren't working and you need to check for missing data at point locations
#rasValue <- extract(predictors_Current, spData_pro)
#write.csv(rasValue, file="rasvalue.csv")

plotCor(data.SWD , method="spearman", cor_th=0.7)


bg_var_sel <- prepareSWD(species=unique(spData_pro$SNAME), a=coords_bg_thin, env=predictors_Current)
plotCor(bg_var_sel, method="spearman", cor_th=0.7)
corVar(bg_var_sel, method="spearman", cor_th=0.7)

#################################################################################################
# run the models ################################################################################

# create folds for crossvalidated model
folds <- randomFolds(data.SWD, k=4, only_presence=TRUE) #, seed=5

for(i in 1:length(ModelMethods)){
  cat("--------------------------------------------------------\n")
  cat(paste("Running a ", ModelMethods[i], " model for ", unique(spData_pro$SNAME),".\n", sep=""))

  cv_model <- train(method=ModelMethods[i], data=data.SWD, folds=folds)  # model <- train(method = "Maxent", data = data, fc = "lh", reg = 0.5, iter = 700)
  cv_model
  
  # calculating 
  cat("- variable importance\n")
  if(ModelMethods[i]=="Maxent"){ 
    md_tssPre <- tss(cv_model)
    md_aucPre <- auc(cv_model)
    cat("Testing TSS before: ", md_tssPre, "\n")
    #vs <- varSel(cv_model, metric="auc", test=NULL, bg4cor=bg_var_sel, method="spearman", cor_th=0.7, permut=10, use_pc=TRUE) # NOTE: This may take a long time...  long time...
    vs <- reduceVar(cv_model, 5, metric="auc", use_pc=TRUE)
    cat("The following variables were used in the final model:", names(vs@data@data), "\n")
    md_tssPost <- tss(vs, test=TRUE)
    md_aucPost <- auc(vs, test=TRUE)
    threshold1 <- as.data.frame(thresholds(vs@models[[1]]))
    threshold2 <- as.data.frame(thresholds(vs@models[[2]]))
    threshold3 <- as.data.frame(thresholds(vs@models[[3]]))
    threshold4 <- as.data.frame(thresholds(vs@models[[4]]))
    thresholds <- bind_rows(threshold1,threshold2,threshold3,threshold4)
    names(thresholds)[2] <- "value" #there is a space in this column name for some reason
    thresholds <- thresholds[thresholds$Threshold=="Minimum training presence",]
    Thresholdmean_minTrainPres <- mean(thresholds$value)
    Thresholdsd_minTrainPres <- sd(thresholds$value)
    cat("Testing TSS after: ", md_tssPost, "\n")    
    # calculate variable importance, this is written to the sqlite db later in the script
    vi <- maxentVarImp(vs)
    plotVarImp(vi[, 1:2])
  } else if(ModelMethods[i]=="RF"|ModelMethods[i]=="BRT") {
    # jacknifing approach to variable selection
    md_tssPre <- tss(cv_model)
    md_aucPre <- auc(cv_model)
    cat("Testing TSS before: ", md_tssPre, "\n")
    vs <- varSel(cv_model, metric="auc", test=NULL, bg4cor=bg_var_sel, method="spearman", cor_th=0.7, permut=10)
    cat("The following variables were used in the final model:", names(vs@data@data), "\n")
    md_tssPost <- tss(vs, test=TRUE)
    md_aucPost <- auc(vs, test=TRUE)
    threshold1 <- as.data.frame(thresholds(vs@models[[1]]))
    threshold2 <- as.data.frame(thresholds(vs@models[[2]]))
    threshold3 <- as.data.frame(thresholds(vs@models[[3]]))
    threshold4 <- as.data.frame(thresholds(vs@models[[4]]))
    thresholds <- bind_rows(threshold1,threshold2,threshold3,threshold4)
    names(thresholds)[2] <- "value" #there is a space in this column name for some reason
    thresholds <- thresholds[thresholds$Threshold=="Minimum training presence",]
    Thresholdmean_minTrainPres <- mean(thresholds$value)
    Thresholdsd_minTrainPres <- sd(thresholds$value)
    cat("Testing TSS after: ", md_tssPost, "\n")
    # calculate variable importance, this is written to the sqlite db later in the script
    vi <- varImp(vs)
    plotVarImp(vi[, 1:2])
  } else {
    cat("No valid variable importance method exists...")
  }

  # insert some model run metadata
  cat("Inserting metadata into the database")
  sf_metadata <- data.frame("sp_code"=sp_code, "model_run_name"=model_run_name, "model_type"=ModelMethods[i], "modeller"=modeller, "TrainingPoints"=md_ptTraining, "TrainingPoints_thinned"=md_ptTrainingThinned, "BackgroundPoints"=md_bg, "BackgroundPoints_thinned"=md_bgThinned, "AUCpre"=md_aucPre, "AUCpost"=md_aucPost, "TSSpre"=md_tssPre, "TSSpost"=md_tssPost, "Thresholdmean_minTrainPres"=Thresholdmean_minTrainPres, "Thresholdsd_minTrainPres"=Thresholdsd_minTrainPres)
  db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
  SQLquery <- paste("INSERT INTO model_runs (", paste(names(sf_metadata), collapse = ',') ,") VALUES (",paste(sQuote(sf_metadata[1,]), collapse = ','),");", sep="") 
  dbExecute(db_cem, SQLquery )
  dbDisconnect(db_cem)

  # insert the variable importance into the database
  md_vi <- data.frame(Variable=character(),
                      Percent_contribution=double(),
                      Permutation_importance=double(),
                      sd=double(),
                      stringsAsFactors=FALSE)
  md_vi <- bind_rows(md_vi, vi)
  md_vi <- cbind("model_run_name"=model_run_name, "model_type"=ModelMethods[i], md_vi)
  for(h in 1:nrow(md_vi)){
    db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
    SQLquery <- paste("INSERT INTO ImpVar (", paste(names(md_vi), collapse = ',') ,") VALUES (",paste(sQuote(md_vi[h,]), collapse = ','),");", sep="") 
    dbExecute(db_cem, SQLquery )
    dbDisconnect(db_cem)
  }
  
  # predict the model to the current env predictors
  cat("- predicting the model to the current env predictors\n")
  timeframe <- "current"
  if(ModelMethods[i]=="Maxent"){
    map <- predict(vs, data=predictors_Current, type="cloglog")
  } else if(ModelMethods[i]=="RF"|ModelMethods[i]=="BRT") {
    map <- predict(vs, data=predictors_Current)
  } else {
    cat("No valid model predictor method exists...\n")
  }
  plotPred(map) # plot and write the current map
  rasnameCurrent <- here::here("_data","species",sp_code,"output",paste(model_run_name, "_", ModelMethods[i], "_", timeframe, ".tif", sep=""))
  writeRaster(map, rasnameCurrent, "GTiff", overwrite=TRUE)
  #predict_current_fn
  db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
  SQLquery <- paste("UPDATE model_runs SET predict_current_fn = ", sQuote(rasnameCurrent), " WHERE model_run_name = ", sQuote(model_run_name), " AND model_type = ", sQuote(ModelMethods[i]), sep="") 
  dbSendStatement(db_cem, SQLquery)
  dbDisconnect(db_cem)
  
  
  # predict the model to the future env predictors
  cat("- predicting the model to the future env predictors\n")
  timeframe <- "future"
  if(ModelMethods[i]=="Maxent"){
    map <- predict(vs, data=predictors_Future4.5, type="cloglog")
  } else if(ModelMethods[i]=="RF"|ModelMethods[i]=="BRT") {
    map <- predict(vs, data=predictors_Future4.5)
  } else {
    cat("No valid model predictor method exists...\n")
  }
  plotPred(map) # plot and write the future map
  rasnameFuture <- here::here("_data","species",sp_code,"output",paste(model_run_name, "_", ModelMethods[i], "_", timeframe, ".tif", sep=""))
  writeRaster(map, rasnameFuture, "GTiff", overwrite=TRUE)
  # insert prediction file names into the database
  
  cat("Inserting more metadata into the database\n")
  
  #predict_future_fn
  db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
  SQLquery <- paste("UPDATE model_runs SET predict_future_fn = ", sQuote(rasnameFuture), " WHERE model_run_name = ", sQuote(model_run_name), " AND model_type = ", sQuote(ModelMethods[i]), sep="") 
  dbSendStatement(db_cem, SQLquery)
  dbDisconnect(db_cem)
  
  
  # predict the model to the future env predictors--using 8.5 scenario
  cat("- predicting the model to the future 8.5 env predictors\n")
  timeframe <- "future8.5"
  if(ModelMethods[i]=="Maxent"){
    map8.5 <- predict(vs, data=predictors_Future8.5, type="cloglog")
  } else if(ModelMethods[i]=="RF"|ModelMethods[i]=="BRT") {
    map8.5 <- predict(vs, data=predictors_Future8.5)
  } else {
    cat("No valid model predictor method exists...\n")
  }
  plotPred(map8.5) # plot and write the future map
  rasnameFuture <- here::here("_data","species",sp_code,"output",paste(model_run_name, "_", ModelMethods[i], "_", timeframe, ".tif", sep=""))
  writeRaster(map, rasnameFuture, "GTiff", overwrite=TRUE)
  
  # insert prediction file names into the database
  cat("Inserting more metadata into the database\n")
  
  #predict_future_fn
  db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
  SQLquery <- paste("UPDATE model_runs SET predict_future_fn85 = ", sQuote(rasnameFuture), " WHERE model_run_name = ", sQuote(model_run_name), " AND model_type = ", sQuote(ModelMethods[i]), sep="") 
  dbSendStatement(db_cem, SQLquery)
  dbDisconnect(db_cem)
  
  # cleanup
  rm(cv_model, vs, vi)
  cat(paste("Finished with the ", ModelMethods[i], " model for ", unique(spData_pro$SNAME),".\n", sep=""))
}

