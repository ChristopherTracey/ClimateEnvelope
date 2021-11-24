######################################
#Build stacked ensemble model ########
#################################################################

#select whether you are running the ensembling within the same session as the individual modeling occurred, or whether are you doing the ensembles in a separate session and starting with a blank workspace

# enter "start fresh" here if you need to select the species and the model run by hand
runtype <- "start fresh" #"start fresh" 

#if you are starting this section separately from running a set of models, then you can manually define your model_run_name here
if (runtype=="start fresh") {
  
sp_code <- "lupipere" ######## manually set focal species HERE
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
model_runs <- dbGetQuery(db_cem, paste("SELECT model_run_name FROM MODEL_RUNS WHERE sp_code = " , sQuote(sp_code), sep="") )
model_runs <- model_runs[,1] #switching from a dataframe to a vector
model_runs
model_run_name <- model_runs[23] ######### manually select the model run name you want to use for ensembling HERE
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
SQLquery <- paste("SELECT model_run_name, model_type, predict_current_fn, predict_future_fn, Thresholdmean_minTrainPres, Thresholdsd_minTrainPres, TSSpost FROM MODEL_RUNS WHERE model_run_name = ", sQuote(model_run_name)) 
model_metadata <- dbGetQuery(db_cem, SQLquery)
dbDisconnect(db_cem)
model_metadata <- unique(model_metadata)
}  else {  

# get model output names from metadata
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
SQLquery <- paste("SELECT model_run_name, model_type, predict_current_fn, predict_future_fn, Thresholdmean_minTrainPres, Thresholdsd_minTrainPres, TSSpost FROM MODEL_RUNS WHERE model_run_name = ", sQuote(model_run_name)) 
model_metadata <- dbGetQuery(db_cem, SQLquery)
dbDisconnect(db_cem)
model_metadata <- unique(model_metadata) # git it down to one row as I have it write out two rows for some reason..
  }

#####
#binarize individual rasters, based on the minimum training presence threshold

# load the current and future rasters
BRT_current <- raster(model_metadata[which(model_metadata$model_type=="BRT"),"predict_current_fn"])
Maxent_current <- raster(model_metadata[which(model_metadata$model_type=="Maxent"),"predict_current_fn"])
RF_current <- raster(model_metadata[which(model_metadata$model_type=="RF"),"predict_current_fn"])

BRT_future <- raster(model_metadata[which(model_metadata$model_type=="BRT"),"predict_future_fn"])
Maxent_future <- raster(model_metadata[which(model_metadata$model_type=="Maxent"),"predict_future_fn"])
RF_future <- raster(model_metadata[which(model_metadata$model_type=="RF"),"predict_future_fn"])

#threshold values
Maxent_t <- model_metadata[which(model_metadata$model_type=="Maxent"),"Thresholdmean_minTrainPres"] 
BRT_t <- model_metadata[which(model_metadata$model_type=="BRT"),"Thresholdmean_minTrainPres"]
RF_t <- model_metadata[which(model_metadata$model_type=="RF"),"Thresholdmean_minTrainPres"]

# function to binerize the MaxEnt model
bin_M <- function(x) {
  ifelse(x <= Maxent_t, 0,
         ifelse(x >  Maxent_t, 1, NA)) }
# function to binerize the BRT model
bin_BRT <- function(x) {
  ifelse(x <= BRT_t, 0,
         ifelse(x >  BRT_t, 1, NA)) }
# function to binerize the RF model
bin_RF <- function(x) {
  ifelse(x <= RF_t, 0,
         ifelse(x >  RF_t, 1, NA)) }

# binarize the current rasters
Maxent_current_bin <- calc(Maxent_current, fun=bin_M)
BRT_current_bin <- calc(BRT_current, fun=bin_BRT)
RF_current_bin <- calc(RF_current, fun=bin_RF)
current_bin <- stack(Maxent_current_bin, BRT_current_bin, RF_current_bin)
current_bin_s <- calc(current_bin, sum)

# stack and average the current maps, weighting by TSS
current <- stack(BRT_current, Maxent_current, RF_current)
current_wm <- weighted.mean(current, w=model_metadata$TSSpost) #weighted mean of the three current models, w/ TSS used to weight

# binarize the future rasters
Maxent_future_bin <- calc(Maxent_future, fun=bin_M)
BRT_future_bin <- calc(BRT_future, fun=bin_BRT)
RF_future_bin <- calc(RF_future, fun=bin_RF)
future_bin <- stack(Maxent_future_bin, BRT_future_bin, RF_future_bin)
future_bin_s <- calc(future_bin, sum)

###########################################
#### Three potential ensemble outputs #####
###########################################

# 1. Continous, stacked, weighted, and averaged 

# stack and average the future maps, weighting by TSS
future <- stack(BRT_future, Maxent_future, RF_future)
future_wm <- weighted.mean(future, w=model_metadata$TSSpost) #weighted mean of the three current models, w/ TSS used to weight. This is a continuous model

plot(future_wm)

# 2. Majority consensus model (sum of 3 binary models, retaining only cells where 2 of 3 models predict presence)
#re-binning the binary consensus model

bin_fut <- function(x) {
  ifelse(x <=  1, 0,
         ifelse(x >  1, 1, NA)) }
future_bin2_s <- calc(future_bin_s, fun=bin_fut)

plot(future_bin2_s) #just the majority consensus points (2 of 3 models agree)

#3. Averaged binned ensemble model (use mean of threshold values to threshold the mean ensemble future model)

mean_threshold <- mean(c(Maxent_t, BRT_t, RF_t))

bin_fut2 <- function(x) {
  ifelse(x <=  mean_threshold, 0,
         ifelse(x >  mean_threshold, 1, NA)) }

future_m <- calc(future, fun=mean) #unweighted mean of three models
current_m <- calc(current, fun=mean) #unweighted mean of three models

future_bin3_s <- calc(future_m, fun=bin_fut2) #this is the thresholded version of the three models averaged, and then the threshold is set by the average. I think this is a thing you can do?
current_bin3_s <- calc(current_m, fun=bin_fut2) #for the current model

plot(current_bin3_s)
plot(future_bin3_s)

##### Save selected ensemble model approach #####
# which model to save and write out as .tif (default is ensemble approach #3)

# create an output directory if it doesn't exist
ifelse(!dir.exists(here::here("_data","thresholded_ensemble_models")), dir.create(here::here("_data","thresholded_ensemble_models")), FALSE)
ensemble_path <- here::here("_data","thresholded_ensemble_models",paste(model_run_name, "_ensemble",".tif", sep=""))

output_model <- future_bin3_s #set which ensemble model you want to export (default is the third option)

writeRaster(output_model, ensemble_path, overwrite=TRUE)

# build CEM expand - contract - stable map

future_bin3_s <- reclassify(future_bin3_s, c(-Inf, .25, 0, .25, 2, 2)) #reclassify the future to a 0,2 raster

cem_cf <- stack(current_bin3_s, future_bin3_s)
cem_cv_s <- calc(cem_cf, fun=sum)

plot(cem_cv_s, legend = FALSE, col = rev(terrain.colors(4)))
legend("topright", legend = c("Null", "Contracting", "Expanding", "Stable"), fill = rev(terrain.colors(4)))

