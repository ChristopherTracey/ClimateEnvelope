library(here)
library(rasterVis)
library(rgdal)
library(grid)
library(scales)
library(viridis)
library(ggthemes) # theme_map()

options(useFancyQuotes=FALSE) # needed to make sure SQL queries work as well as they could

source(here::here("0_runCEM.r"))

######################################
#Build stacked ensemble model ########
#################################################################

sp_code <- sp_code
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
model_runs <- dbGetQuery(db_cem, paste("SELECT model_run_name FROM MODEL_RUNS WHERE sp_code = " , sQuote(sp_code), sep="") )
dbDisconnect(db_cem)
model_runs <- unique(model_runs[,1]) #switching from a dataframe to a vector

if(length(model_runs)==1){
  model_run <- model_runs
} else if(length(model_runs)>1){
  cat("Multiple models found, please select the number of the model you want to use:\n")
  cat(model_runs)
  n <- 2
  model_run <- model_runs[n]
} else {
  cat("No model found. Was it run yet?")
}

db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
SQLquery <- paste("SELECT model_run_name, model_type, predict_current_fn, predict_future_fn, predict_future_fn, predict_future_fn85, Thresholdmean_minTrainPres, Thresholdsd_minTrainPres, TSSpost FROM MODEL_RUNS WHERE model_run_name = ", sQuote(model_run)) 
model_metadata <- dbGetQuery(db_cem, SQLquery)
dbDisconnect(db_cem)
model_metadata <- unique(model_metadata)

model_run_name <- model_run

###### binarize individual rasters, based on the minimum training presence threshold
# load the current rasters
BRT_current <- raster(model_metadata[which(model_metadata$model_type=="BRT"),"predict_current_fn"])
Maxent_current <- raster(model_metadata[which(model_metadata$model_type=="Maxent"),"predict_current_fn"])
RF_current <- raster(model_metadata[which(model_metadata$model_type=="RF"),"predict_current_fn"])
# load the future 4.5 rasters
BRT_future45 <- raster(model_metadata[which(model_metadata$model_type=="BRT"),"predict_future_fn"])
Maxent_future45 <- raster(model_metadata[which(model_metadata$model_type=="Maxent"),"predict_future_fn"])
RF_future45 <- raster(model_metadata[which(model_metadata$model_type=="RF"),"predict_future_fn"])
# load the future 8.5 rasters
BRT_future85 <- raster(model_metadata[which(model_metadata$model_type=="BRT"),"predict_future_fn85"])
Maxent_future85 <- raster(model_metadata[which(model_metadata$model_type=="Maxent"),"predict_future_fn85"])
RF_future85 <- raster(model_metadata[which(model_metadata$model_type=="RF"),"predict_future_fn85"])

#threshold values
Maxent_t <- model_metadata[which(model_metadata$model_type=="Maxent"),"Thresholdmean_minTrainPres"]
BRT_t <- model_metadata[which(model_metadata$model_type=="BRT"),"Thresholdmean_minTrainPres"]
RF_t <- model_metadata[which(model_metadata$model_type=="RF"),"Thresholdmean_minTrainPres"]

#state outlines
states <- arc.open("https://maps.waterlandlife.org/arcgis/rest/services/BaseLayers/Boundaries/FeatureServer/3")
states <- arc.select(states)
states <- arc.data2sf(states)
states <- states[states$NAME %in% c("Pennsylvania", "New York", "Delaware", "West Virginia", "Virginia", "Ohio", "Maryland", "New Jersey"),]
states <- st_transform(states, crs=st_crs(BRT_current))
states <- states %>% group_by(NAME) %>% summarize() # dissolve singleparts to multipart (e.g. long island)

#Species points
Sp_path <- here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.shp"))
spData <- arc.open(Sp_path)
spData <- arc.select(spData)
spData <- arc.data2sf(spData)
spData_pro <- st_transform(spData, crs=crs(BRT_current))

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
ggsave(filename=here::here("_data", "species", sp_code, "output", paste(model_run_name, "_current_wm_pts",".png", sep="")), plot=current_wm_pts, device='png', dpi=300)

#create map and export WITHOUT spp points overlaid
current_wm <- ggplot() + geom_raster(data=current_wmdf, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_sf(data=states, fill=NA, color="black") +
  scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
ggsave(filename=here::here("_data", "species", sp_code, "output", paste(model_run_name, "_current_wm",".png", sep="")), plot=current_wm, device='png', dpi=300)

# binarize the future rasters (4.5, the baseline reporting scenario we're using)
Maxent_future45_bin <- calc(Maxent_future45, fun=bin_M)
BRT_future45_bin <- calc(BRT_future45, fun=bin_BRT)
RF_future45_bin <- calc(RF_future45, fun=bin_RF)
future45_bin <- stack(Maxent_future45_bin, BRT_future45_bin, RF_future45_bin)
future45_bin_s <- calc(future45_bin, sum)

#binarize the future 8.5 rasters
Maxent_future_bin85 <- calc(Maxent_future85, fun=bin_M)
BRT_future_bin85 <- calc(BRT_future85, fun=bin_BRT)
RF_future_bin85 <- calc(RF_future85, fun=bin_RF)
future_bin85 <- stack(Maxent_future_bin85, BRT_future_bin85, RF_future_bin85)
future_bin_s85 <- calc(future_bin85, sum)

# stack and average the future maps, weighting by TSS
future45 <- stack(BRT_future45, Maxent_future45, RF_future45)
future45_wm <- weighted.mean(future45, w=model_metadata$TSSpost) #weighted mean of the three current models, w/ TSS used to weight. This is a continuous model
future45_wmdf <- as.data.frame(future45_wm, xy = TRUE) %>% na.omit()
names(future45_wmdf) <- names_rdf

# stack and average the future maps for the 8.5 scenario, weighting by TSS
future85 <- stack(BRT_future85, Maxent_future85, RF_future85)
future85_wm <- weighted.mean(future85, w=model_metadata$TSSpost) #weighted mean of the three current models, w/ TSS used to weight. This is a continuous model
future85_wmdf <- as.data.frame(future85_wm, xy = TRUE) %>% na.omit()
names(future85_wmdf) <- names_rdf

plot(future85_wm)

#create map and export WITH spp points overlaid
future45_wm_pts <- ggplot() + 
  geom_raster(data=future45_wmdf, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_point(data=sp_pts, aes(x=X, y=Y),shape=3) + geom_sf(data=states, fill=NA, color="black") +
  scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
ggsave(filename=here::here("_data", "species", sp_code, "output", paste(model_run_name, "_future_45_wm_pts",".png", sep="")), plot=future45_wm_pts, device='png', dpi=300)

#create map and export WITHOUT spp points overlaid
future45_wm <- ggplot() + 
  geom_raster(data=future45_wmdf, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_sf(data=states, fill=NA, color="black") +
  scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
ggsave(filename=here::here("_data", "species", sp_code, "output", paste(model_run_name, "_future_45_wm",".png", sep="")), plot=future45_wm, device='png', dpi=300)

#and do the same for the 8.5 future scenario
future85_wm_pts <- ggplot() + 
  geom_raster(data=future85_wmdf, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_point(data=sp_pts, aes(x=X, y=Y),shape=3) + geom_sf(data=states, fill=NA, color="black") +
  scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
ggsave(filename=here::here("_data", "species", sp_code, "output", paste(model_run_name, "_future_85_wm_pts",".png", sep="")), plot=future85_wm_pts, device='png', dpi=300)

#create map and export WITHOUT spp points overlaid
future85_wm <- ggplot() + 
  geom_raster(data=future85_wmdf, aes(x=x, y=y, fill=Likelihood), alpha=0.7) +
  geom_sf(data=states, fill=NA, color="black") +
  scale_fill_viridis_c(limits = c(0, 1)) + theme_void() + theme(legend.position = "bottom")
ggsave(filename=here::here("_data", "species", sp_code, "output", paste(model_run_name, "_future_85_wm",".png", sep="")), plot=future85_wm, device='png', dpi=300)

#Averaged binned ensemble model (use mean of threshold values to threshold the mean ensemble future model)
mean_threshold <- mean(c(Maxent_t, BRT_t, RF_t))

bin_fut2 <- function(x) {
  ifelse(x <=  mean_threshold, 0,
         ifelse(x >  mean_threshold, 1, NA)) }

current_m <- calc(current, fun=mean) # unweighted mean of three models
future45_m <- calc(future45, fun=mean) # unweighted mean of three models
future85_m <- calc(future85, fun=mean) # unweighted mean of three models

# this is the thresholded version of the three models averaged, and then the threshold is set by the average. I think this is a thing you can do?
current_bin3_s <- calc(current_m, fun=bin_fut2) 
future45_bin3_s <- calc(future45_m, fun=bin_fut2) 
future85_bin3_s <- calc(future85_m, fun=bin_fut2)

plot(current_bin3_s)
plot(future45_bin3_s)
plot(future85_bin3_s)

##### Save ensemble binned models in species output folder#####
writeRaster(current_bin3_s, here::here("_data", "species", sp_code, "output", paste(model_run_name, "_ensemble_thresh", "_", "current",".tif", sep="")), overwrite=TRUE)
writeRaster(future45_bin3_s, here::here("_data", "species", sp_code, "output", paste(model_run_name, "_ensemble_thresh", "_", "future45",".tif", sep="")), overwrite=TRUE)
writeRaster(future85_bin3_s, here::here("_data", "species", sp_code, "output", paste(model_run_name, "_ensemble_thresh", "_", "future85",".tif", sep="")), overwrite=TRUE)

# build CEM expand - contract - stable map, using the 4.5 scenario
future45_bin3_s <- reclassify(future45_bin3_s, c(-Inf, .25, 0, .25, 2, 2)) #reclassify the future to a 0,2 raster

cem_cf <- stack(current_bin3_s, future45_bin3_s)
cem_cv_s <- calc(cem_cf, fun=sum)
cem_cv_sdf <- as.data.frame(cem_cv_s, xy = TRUE) %>% na.omit()
names(cem_cv_sdf) <- c("x", "y", "Range.Shift")
cem_cv_sdf$Range.Shift <- factor(cem_cv_sdf$Range.Shift, levels=c(0,1,2,3), labels=c("Null", "Contracting", "Expanding", "Stable" ))
rm(cem_cf)

# write the cem to a tif
writeRaster(cem_cv_s, here::here("_data", "species", sp_code, "output", paste(model_run_name, "_cem",".tif", sep="")), overwrite=TRUE)

# summarize by area
fr <- rasterize(states, cem_cv_s) # rasterize the states to mask outside project area
lr <- mask(x=cem_cv_s, mask=fr) # mask by the project area
rm(fr)

tabFunc<-function(indx, extracted, region, regname) {
  dat <- as.data.frame(table(extracted[[indx]]))
  dat$name <- region[[regname]][[indx]]
  return(dat)
}

ext <- extract(lr, states, method='simple')
cem_areaTab <- lapply(seq(ext), tabFunc, ext, states, "NAME")
cem_areaTab <- do.call("rbind", cem_areaTab) # convert to df
cem_areaTab$Var1 <- factor(cem_areaTab$Var1, levels=c(0,1,2,3), labels=c("Null", "Contracting", "Expanding", "Stable" ))
cem_areaTab_wide <- cem_areaTab %>%
  group_by(name) %>% # group by region
  mutate(totcells=sum(Freq), # how many cells overall
         percent.area=round(100*Freq/totcells,2)) %>% #cells by shift/total cells
  dplyr::select(-c(Freq, totcells)) %>% # there is a select func in raster so need to specify
  spread(key=Var1, value=percent.area, fill=0)

# make the plots for the metadata
custom_fill_pal <- c(Contracting="#ECB176", Stable="#00A600", Expanding="#E6E600", Null="#F2F2F2")
#THRESHOLD CONTRACT EXPAND STABLE MAP, NO POINTS
threshold_plot <- ggplot() + 
  geom_raster(data=cem_cv_sdf, aes(x=x, y=y, fill=Range.Shift)) + 
  scale_fill_manual(values=custom_fill_pal) + 
  geom_sf(data=states, fill=NA, color="black") + 
  theme_void()
ggsave(filename=here::here("_data", "species", sp_code, "output", paste(model_run_name, "_thresholdPlot",".png", sep="")), plot=threshold_plot, device='png', dpi=300)

#THRESHOLD CONTRACT EXPAND STABLE MAP, WITH SPECIES POINTS
threshold_plot_SP <- ggplot() + 
  geom_raster(data=cem_cv_sdf, aes(x=x, y=y, fill=Range.Shift)) + 
  scale_fill_manual(values= custom_fill_pal) + 
  geom_point(data=sp_pts, aes(x=X, y=Y), shape=3) + 
  geom_sf(data=states, fill=NA, color="black") + 
  theme_void()
ggsave(filename=here::here("_data", "species", sp_code, "output", paste(model_run_name, "_thresholdPlot_pts",".png", sep="")), plot=threshold_plot_SP, device='png', dpi=300)





# 
# #### Notes, other stuff, doesn't need to run ####
# 
# #### alternative ensemble approaches (not currently using)
# 
# # 2. Majority consensus model (sum of 3 binary models, retaining only cells where 2 of 3 models predict presence)
# #re-binning the binary consensus model
# 
# bin_fut <- function(x) {
#   ifelse(x <=  1, 0,
#          ifelse(x >  1, 1, NA)) }
# future_bin2_s <- calc(future_bin_s, fun=bin_fut)
# 
# plot(future_bin2_s) #just the majority consensus points (2 of 3 models agree)
