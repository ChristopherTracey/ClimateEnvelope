library(here)
library(arcgisbinding)
arc.check_product()
library(ggplot2)
library(RSQLite)
library(sf)
library(tidyverse)
library(SDMtune)
library(reshape2)
library(raster)
library(remotes)
library(dplyr)
library(lubridate)
library(rasterVis)
library(rgdal)
library(grid)
library(scales)
library(viridis)
library(ggthemes) # theme_map()

options(useFancyQuotes=FALSE) # needed to make sure SQL queries work as well as they could

######################################
#Build stacked ensemble model ########
#################################################################

#select whether you are running the ensembling within the same session as the individual modeling occurred, or whether are you doing the ensembles in a separate session and starting with a blank workspace

# enter "start fresh" here if you need to select the species and the model run by hand
runtype <- "continue" #"start fresh" #continue

#file path to save .tiff and geotiff outputs from ensemble modeling
map_path <- here::here("_data", "species", sp_code, "output")

#if you are starting this section separately from running a set of models, then you can manually define your model_run_name here
if (runtype=="start fresh") {

model_run_name <- "abiebals_20211222_135814" #input the model run name you need
sp_code <- "abiebals" # which species are you working with? Abies balsamifera  
  
sp_code <- sp_code
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
model_runs <- dbGetQuery(db_cem, paste("SELECT model_run_name FROM MODEL_RUNS WHERE sp_code = " , sQuote(sp_code), sep="") )
model_runs <- model_runs[,1] #switching from a dataframe to a vector
model_runs
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
SQLquery <- paste("SELECT model_run_name, model_type, predict_current_fn, predict_future_fn, predict_future_fn, predict_future_fn85, Thresholdmean_minTrainPres, Thresholdsd_minTrainPres, TSSpost FROM MODEL_RUNS WHERE model_run_name = ", sQuote(model_run_name)) 
model_metadata <- dbGetQuery(db_cem, SQLquery)
dbDisconnect(db_cem)
model_metadata <- unique(model_metadata)
}  else {  

# get model output names from metadata
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
SQLquery <- paste("SELECT model_run_name, model_type, predict_current_fn, predict_future_fn, predict_future_fn85, Thresholdmean_minTrainPres, Thresholdsd_minTrainPres, TSSpost FROM MODEL_RUNS WHERE model_run_name = ", sQuote(model_run_name)) 
model_metadata <- dbGetQuery(db_cem, SQLquery)
dbDisconnect(db_cem)
model_metadata <- unique(model_metadata) # git it down to one row as I have it write out two rows for some reason..
  }

#####
crs_map <- "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#binarize individual rasters, based on the minimum training presence threshold

# load the current and future rasters
BRT_current <- raster(model_metadata[which(model_metadata$model_type=="BRT"),"predict_current_fn"])
crs(BRT_current) <- crs_map
Maxent_current <- raster(model_metadata[which(model_metadata$model_type=="Maxent"),"predict_current_fn"])
crs(Maxent_current) <- crs_map
RF_current <- raster(model_metadata[which(model_metadata$model_type=="RF"),"predict_current_fn"])
crs(RF_current) <- crs_map

BRT_future <- raster(model_metadata[which(model_metadata$model_type=="BRT"),"predict_future_fn"])
crs(BRT_future) <- crs_map
Maxent_future <- raster(model_metadata[which(model_metadata$model_type=="Maxent"),"predict_future_fn"])
crs(Maxent_future) <- crs_map
RF_future <- raster(model_metadata[which(model_metadata$model_type=="RF"),"predict_future_fn"])
crs(RF_future) <- crs_map

BRT_future85 <- raster(model_metadata[which(model_metadata$model_type=="BRT"),"predict_future_fn85"])
crs(BRT_future85) <- crs_map
Maxent_future85 <- raster(model_metadata[which(model_metadata$model_type=="Maxent"),"predict_future_fn85"])
crs(Maxent_future85) <- crs_map
RF_future85 <- raster(model_metadata[which(model_metadata$model_type=="RF"),"predict_future_fn85"])
crs(RF_future85) <- crs_map

#threshold values
Maxent_t <- model_metadata[which(model_metadata$model_type=="Maxent"),"Thresholdmean_minTrainPres"]
BRT_t <- model_metadata[which(model_metadata$model_type=="BRT"),"Thresholdmean_minTrainPres"]
RF_t <- model_metadata[which(model_metadata$model_type=="RF"),"Thresholdmean_minTrainPres"]

#state outlines
states <- arc.open("https://maps.waterlandlife.org/arcgis/rest/services/BaseLayers/Boundaries/FeatureServer/3")
states <- arc.select(states)
states <- arc.data2sf(states)
states <- states[states$NAME %in% c("Pennsylvania","New York", "Delaware", "West Virginia", "Virginia", "Ohio", "Maryland", "New Jersey"),]
states <- st_transform(states, crs=st_crs(BRT_current))

#Species points
Sp_path <- here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.shp"))
spData <- arc.open(Sp_path)
spData <- arc.select(spData)
spData <- arc.data2sf(spData)
spData_pro <- st_transform(spData, crs=crs(crs_map))

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

plot(current_wm)
current_wmdf <- as.data.frame(current_wm, xy = TRUE) %>% na.omit()
names_rdf <- c("x","y","Likelihood")
names(current_wmdf) <- names_rdf
sp_pts <- as.data.frame(st_coordinates(spData_pro))

#create map and export WITH spp points overlaid
current_wm_pts <- ggplot() + geom_raster(data=current_wmdf, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_point(data=sp_pts, aes(x=X, y=Y),shape=3) + geom_sf(data=states, fill=NA, color="black") +
scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
current_wm_pts

ggsave(filename="current_wm_pts.png", plot=last_plot(), path = map_path, device='png', dpi=300)

#create map and export WITHOUT spp points overlaid
current_wm <- ggplot() + geom_raster(data=current_wmdf, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_sf(data=states, fill=NA, color="black") +
  scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
current_wm

ggsave(filename="current_wm.png", plot=last_plot(), path = map_path, device='png', dpi=300)

# binarize the future rasters (4.5, the baseline reporting scenario we're using)
Maxent_future_bin <- calc(Maxent_future, fun=bin_M)
BRT_future_bin <- calc(BRT_future, fun=bin_BRT)
RF_future_bin <- calc(RF_future, fun=bin_RF)
future_bin <- stack(Maxent_future_bin, BRT_future_bin, RF_future_bin)
future_bin_s <- calc(future_bin, sum)

#binarize the future 8.5 rasters
Maxent_future_bin85 <- calc(Maxent_future85, fun=bin_M)
BRT_future_bin85 <- calc(BRT_future85, fun=bin_BRT)
RF_future_bin85 <- calc(RF_future85, fun=bin_RF)
future_bin85 <- stack(Maxent_future_bin85, BRT_future_bin85, RF_future_bin85)
future_bin_s85 <- calc(future_bin85, sum)

# stack and average the future maps, weighting by TSS
future <- stack(BRT_future, Maxent_future, RF_future)
future_wm <- weighted.mean(future, w=model_metadata$TSSpost) #weighted mean of the three current models, w/ TSS used to weight. This is a continuous model
future_wmdf <- as.data.frame(future_wm, xy = TRUE) %>% na.omit()
names(future_wmdf) <- names_rdf

# stack and average the future maps for the 8.5 scenario, weighting by TSS
future85 <- stack(BRT_future85, Maxent_future85, RF_future85)
future_wm85 <- weighted.mean(future85, w=model_metadata$TSSpost) #weighted mean of the three current models, w/ TSS used to weight. This is a continuous model
future_wmdf85 <- as.data.frame(future_wm85, xy = TRUE) %>% na.omit()
names(future_wmdf85) <- names_rdf

plot(future_wm85)

#create map and export WITH spp points overlaid
future_wm_pts <- ggplot() + geom_raster(data=future_wmdf, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_point(data=sp_pts, aes(x=X, y=Y),shape=3) + geom_sf(data=states, fill=NA, color="black") +
  scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
future_wm_pts

ggsave(filename="future_wm_pts.png", plot=last_plot(), path = map_path, device='png', dpi=300)

#create map and export WITHOUT spp points overlaid
future_wm <- ggplot() + geom_raster(data=future_wmdf, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_sf(data=states, fill=NA, color="black") +
  scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
future_wm

ggsave(filename="future_wm.png", plot=last_plot(), path = map_path, device='png', dpi=300)

#and do the same for the 8.5 future scenario
future_wm_pts85 <- ggplot() + geom_raster(data=future_wmdf85, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_point(data=sp_pts, aes(x=X, y=Y),shape=3) + geom_sf(data=states, fill=NA, color="black") +
  scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
future_wm_pts85

ggsave(filename="future_wm_pts85.png", plot=last_plot(), path = map_path, device='png', dpi=300)

#create map and export WITHOUT spp points overlaid
future_wm85 <- ggplot() + geom_raster(data=future_wmdf85, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_sf(data=states, fill=NA, color="black") +
  scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
future_wm85

ggsave(filename="future_wm85.png", plot=last_plot(), path = map_path, device='png', dpi=300)

#Averaged binned ensemble model (use mean of threshold values to threshold the mean ensemble future model)

mean_threshold <- mean(c(Maxent_t, BRT_t, RF_t))

bin_fut2 <- function(x) {
  ifelse(x <=  mean_threshold, 0,
         ifelse(x >  mean_threshold, 1, NA)) }

future_m <- calc(future, fun=mean) #unweighted mean of three models
current_m <- calc(current, fun=mean) #unweighted mean of three models
future_m85 <- calc(future85, fun=mean)

future_bin3_s <- calc(future_m, fun=bin_fut2) #this is the thresholded version of the three models averaged, and then the threshold is set by the average. I think this is a thing you can do?
current_bin3_s <- calc(current_m, fun=bin_fut2) #for the current model, same deal
future85_bin3_s <- calc(future_m85, fun=bin_fut2)

plot(current_bin3_s)
plot(future_bin3_s)
plot(future85_bin3_s)





##### Save ensemble binned models in species output folder#####
writeRaster(current_bin3_s, map_path, overwrite=TRUE)
writeRaster(future_bin3_s, map_path, overwrite=TRUE)
writeRaster(future85_bin3_s, map_path, overwrite=TRUE)

# build CEM expand - contract - stable map, using the 4.5 scenario

future_bin3_s <- reclassify(future_bin3_s, c(-Inf, .25, 0, .25, 2, 2)) #reclassify the future to a 0,2 raster

cem_cf <- stack(current_bin3_s, future_bin3_s)
cem_cv_s <- calc(cem_cf, fun=sum)

cem_cv_sdf <- as.data.frame(cem_cv_s, xy = TRUE) %>% na.omit()
names_threshdf <- c("x","y","Range.Shift")
names(cem_cv_sdf) <- names_threshdf
cem_cv_sdf$Range.Shift <- as.factor(cem_cv_sdf$Range.Shift)
levels(cem_cv_sdf$Range.Shift) <- c("Null", "Contracting", "Expanding", "Stable")

#THRESHOLD CONTRACT EXPAND STABLE MAP, NO POINTS
custom_fill_pal <- c("#F2F2F2","#ECB176","#E6E600","#00A600")
threshold_plot <- ggplot() + geom_raster(data=cem_cv_sdf, aes(x=x, y=y, fill=Range.Shift)) + scale_fill_manual(values= custom_fill_pal) +  geom_sf(data=states, fill=NA, color="black") + theme_void()
threshold_plot
ggsave(filename="threshold_plot.png", plot=last_plot(), path = map_path, device='png', dpi=300)

#THRESHOLD CONTRACT EXPAND STABLE MAP, WITH SPECIES POINTS
threshold_plot_SP <- ggplot() + geom_raster(data=cem_cv_sdf, aes(x=x, y=y, fill=Range.Shift)) + scale_fill_manual(values= custom_fill_pal) + geom_point(data=sp_pts, aes(x=X, y=Y),shape=3) + geom_sf(data=states, fill=NA, color="black") + theme_void()
threshold_plot_SP

ggsave(filename="threshold_plot_pts.png", plot=last_plot(), path = map_path, device='png', dpi=300)





#### Notes, other stuff, doesn't need to run ####

#### alternative ensemble approaches (not currently using)

# 2. Majority consensus model (sum of 3 binary models, retaining only cells where 2 of 3 models predict presence)
#re-binning the binary consensus model

bin_fut <- function(x) {
  ifelse(x <=  1, 0,
         ifelse(x >  1, 1, NA)) }
future_bin2_s <- calc(future_bin_s, fun=bin_fut)

plot(future_bin2_s) #just the majority consensus points (2 of 3 models agree)
