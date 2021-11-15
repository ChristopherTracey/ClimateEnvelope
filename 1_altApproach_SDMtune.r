
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
st_write(spData, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.shp")), append=FALSE)

# convert species points to a lat/long df  MOVE THIS BELOW???
coords_pres <- data.frame(st_coordinates(spData_pro[,1]))
names(coords_pres) <- c("x","y")
md_ptTraining <- nrow(coords_pres) # get some metadata

# thin out the points to only one per cell
coords_pres_thin <- thinData(coords_pres, predictors_Current)
md_ptTrainingThinned <- nrow(spData) # get some metadata

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
data.SWD <- prepareSWD(species=unique(spData$SNAME), p=coords_pres_thin, a=coords_bg_thin, env=predictors_Current)
data.SWD

plotCor(data.SWD , method="spearman", cor_th=0.7)


bg_var_sel <- prepareSWD(species=unique(spData$SNAME), a=coords_bg_thin, env=predictors_Current)
plotCor(bg_var_sel, method="spearman", cor_th=0.7)
corVar(bg_var_sel, method="spearman", cor_th=0.7)

#################################################################################################
# run the models ################################################################################

# create folds for crossvalidated model
folds <- randomFolds(data.SWD, k=4, only_presence=TRUE) #, seed=5

for(i in 1:length(ModelMethods)){
  cat("--------------------------------------------------------\n")
  cat(paste("Running a ", ModelMethods[i], " model for ", unique(spData$SNAME),".\n", sep=""))

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
  # md_vi <- cbind("model_run_name"=model_run_name, "model_type"=ModelMethods[i], vi)
  # for(h in 1:nrow(md_vi)){
  #   db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
  #   SQLquery <- paste("INSERT INTO ImpVar (", paste(names(md_vi), collapse = ',') ,") VALUES (",paste(sQuote(md_vi[h,]), collapse = ','),");", sep="") 
  #   dbExecute(db_cem, SQLquery )
  #   dbDisconnect(db_cem)
  # }
  
  # predict the model to the current env predictors
  cat("- predicting the model to the current env predictors\n")
  timeframe <- "current"
  if(ModelMethods[i]=="Maxent"){
    map <- predict(vs, data=predictors_Current, type="cloglog")
  } else if(ModelMethods[i]=="RF"|ModelMethods[i]=="BRT") {
    map <- predict(vs, data=predictors_Current)
  } else {
    cat("No valid model predictor method exists...")
  }
  plotPred(map) # plot and write the current map
  rasnameCurrent <- here::here("_data","species",sp_code,"output",paste(model_run_name, "_", ModelMethods[i], "_", timeframe, ".tif", sep=""))
  writeRaster(map, rasnameCurrent, "GTiff", overwrite=TRUE)
  
  # predict the model to the future env predictors
  cat("- predicting the model to the future env predictors\n")
  timeframe <- "future"
  if(ModelMethods[i]=="Maxent"){
    map <- predict(vs, data=predictors_Future, type="cloglog")
  } else if(ModelMethods[i]=="RF"|ModelMethods[i]=="BRT") {
    map <- predict(vs, data=predictors_Future)
  } else {
    cat("No valid model predictor method exists...")
  }
  plotPred(map) # plot and write the future map
  rasnameFuture <- here::here("_data","species",sp_code,"output",paste(model_run_name, "_", ModelMethods[i], "_", timeframe, ".tif", sep=""))
  writeRaster(map, rasnameFuture, "GTiff", overwrite=TRUE)
  
  # insert prediction file names into the database
  cat("Inserting more metadata into the database")
  

  #predict_future_fn
  db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
  SQLquery <- paste("UPDATE model_runs SET predict_current_fn = ", sQuote(rasnameCurrent), " WHERE model_run_name = ", sQuote(model_run_name), " AND model_type = ", sQuote(ModelMethods[i]), sep="") 
  dbSendStatement(db_cem, SQLquery)
#  dbExecute(db_cem, SQLquery )
  dbDisconnect(db_cem)

  # cleanup
  rm(cv_model, vs, vi)
  cat(paste("Finished with the ", ModelMethods[i], " model for ", unique(spData$SNAME),".\n", sep=""))
}

######################################
#Build stacked ensemble model ########
#################################################################

#path to model outputs:
models <- list.files(Model_outputpath) #pull list of files within the model output folder, then use this to select which file to use for which model type

#select which element from models list should be used for each model (change numbers assigned as needed) 
BRT_current <- models[27]
BRT_future <- models[28]
Maxent_current <- models[29]
Maxent_future <- models[30]
RF_current <- models[31]
RF_future <- models[32]

model_v <- c(BRT_current, Maxent_future, RF_current) #just need one row of metadata for a pair of current/future models
model_v2 <- stringr::str_extract(model_v, "[^_]*_[^_]*_[^_]*") #extract everything before the third underscore (to match SQL naming in db)
model_methods <- c("BRT","Maxent","RF")

BRT_currentp <- paste(Model_outputpath,BRT_current,sep="/")
BRT_futurep <- paste(Model_outputpath,BRT_future,sep="/")
Maxent_currentp <- paste(Model_outputpath,Maxent_current,sep="/")
Maxent_futurep <- paste(Model_outputpath,Maxent_future,sep="/")
RF_currentp <- paste(Model_outputpath,RF_current,sep="/")
RF_futurep <- paste(Model_outputpath,RF_future,sep="/")

#pull the associated model metadata
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
SQLquery <- paste("SELECT * FROM MODEL_RUNS WHERE MODEL_RUN_NAME IN (", paste0(sprintf("'%s'", model_v2), collapse = ", "), ")") 
model_metadata <- dbGetQuery(db_cem, SQLquery)
dbDisconnect(db_cem)

#if there are too many rows (i.e. a duplicate model run, fix that mnaually)
#model_metadata <- model_metadata[-1,]
model_metadata <- model_metadata[match(model_methods, model_metadata$model_type),] #sort metadata so it is in the same model type order as the model methods vector

#binarize individual rasters, based on the minimum training presence threshold

#threshold values
Maxent_t <- model_metadata$Thresholdmean_minTrainPres[2]
BRT_t <- model_metadata$Thresholdmean_minTrainPres[1]
RF_t <- model_metadata$Thresholdmean_minTrainPres[3]

bin_M <- function(x) {
  ifelse(x <=  Maxent_t, 0,
         ifelse(x >  Maxent_t, 1, NA)) }

bin_BRT <- function(x) {
  ifelse(x <=  BRT_t, 0,
         ifelse(x >  BRT_t, 1, NA)) }

bin_RF <- function(x) {
  ifelse(x <=  RF_t, 0,
         ifelse(x >  RF_t, 1, NA)) }

#stack and average the current maps, weighting by TSS
BRT <- raster(BRT_currentp)
Maxent <- raster(Maxent_currentp)
RF <- raster(RF_currentp)

#binarize
Maxent_current_bin <- calc(Maxent, fun=bin_M)
BRT_current_bin <- calc(BRT, fun=bin_BRT)
RF_current_bin <- calc(RF, fun=bin_RF)
current_bin <- stack(Maxent_current_bin, BRT_current_bin, RF_current_bin)
current_bin_s <- calc(current_bin, sum)

current <- stack(BRT, Maxent, RF)
current_wm <- weighted.mean(current, w=model_metadata$TSSpost) #weighted mean of the three current models, w/ TSS used to weight

#stack and average the future maps, weighting by TSS
BRT <- raster(BRT_futurep)
Maxent <- raster(Maxent_futurep)
RF <- raster(RF_futurep)

#binarize
Maxent_future_bin <- calc(Maxent, fun=bin_M)
BRT_future_bin <- calc(BRT, fun=bin_BRT)
RF_future_bin <- calc(RF, fun=bin_RF)
future_bin <- stack(Maxent_future_bin, BRT_future_bin, RF_future_bin)
future_bin_s <- calc(future_bin, sum)

future <- stack(BRT, Maxent, RF)
future_wm <- weighted.mean(future, w=model_metadata$TSSpost) #weighted mean of the three current models, w/ TSS used to weight


Maxent_current_bin <- calc(Maxent, fun=bin_M)

#re-binning the binary consensus model
bin_fut <- function(x) {
  ifelse(x <=  2, 0,
         ifelse(x >  2, 1, NA)) }
future_bin2_s <- calc(future_bin_s, fun=bin_fut)
plot(future_bin2_s) #just the full consensus points
