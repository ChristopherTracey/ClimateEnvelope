# GOALS
# - run three models (maxent, rf, and glm)
# - evaluate each individual model
# - create framework to compare and create ensemble models
# - stack models within species groups

library(here)
library(arcgisbinding)
arc.check_product()
arc.check_portal()
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
# library(dismo)
#library(sdm)
#library(rasterVis)

options(useFancyQuotes=FALSE) # needed to make sure SQL queries work as well as they could
source(here::here("helperFunctions.r"))

## STEP 1: Variables and Such ################################
# species code (from lkpSpecies in modelling database. This will be the new folder name containing inputs/outputs)
sp_code <- "utriinte" # Abies balsamifera
# model methods
ModelMethods <- c("Maxent","BRT","RF")
# Modeling database
nm_db_file <- here("_data", "databases", "CEMdata.sqlite")
# species data 
spData_path <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/All_SpPoints/Sp_Points2.shp"
# project area, shapefile or gdb feature class
studyArea <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/Refugia Modeling Boundary/Refugia Modeling Boundary.shp"

# predictor datasets
pathPredictorsCurrent <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_NORM_8110_Bioclim_ASCII/NewVars_masked"
pathPredictorsFuture4.5 <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_ENSEMBLE_rcp45_2050s_Bioclim_ASCII/NewVars_masked"
pathPredictorsFuture8.5 <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_ENSEMBLE_rcp85_2050s_Bioclim_ASCII/NewVars_masked"
# your name
modeller = "Anna Johnson"

#model output path
#Model_outputpath <- here::here(paste("_data/","species","/",sp_code,"/output", sep=""))

# set up projection parameter for use throughout script--albers
#ascproj <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +") #projection of AdaptWest asc layers
#proj <- CRS("+proj=aea +lat_1=40 +lat_2=42 +lat_0=39 +lon_0=-78 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs") #projection for PA


## STEP 2: Run a Model ################################

#source(here::here("1_altApproach_SDMtune.r"))

# # Step 
# 
# source("1_PrepSpeciesData.r")
# source("2_AttributeAndBackground.r")
# source("3a_MaxEnt.r")



# 
# #####
# 
# 
# # loc_scripts is your repository. Make sure your git repository is set to correct branch
# loc_scripts <- here()
# # The main modelling folder for inputs/outputs. All sub-folders are created during the model run (when starting with step 1)
# loc_model <- here("_data", "species")
# 
# # locations file (presence reaches). Provide full path; File is copied to modeling folder and timestamped.
# nm_presFile <- here("_data", "occurrence", paste0(model_species, ".shp"))
# #nm_presFile <- here("_data", "occurrence", paste0(sub("-","_",model_species), ".gpkg"))
# # env vars location [Terrestrial-only variable]
# loc_envVars = here("_data","env_vars","raster", "ras") # "D:/R_tmp/taenmont_20210304_134229"
# # Name of background/envvars sqlite geodatabase, and base table name (2 length vector)
# nm_bkgPts <- c(here("_data","env_vars","tabular", "background_CONUS.sqlite"), "background_pts")
# # HUC spatial data set (shapefile) that is subsetted and used to define modeling area//range
# nm_HUC_file <- here("_data","other_spatial","feature","HUC10.shp")
# # map reference boundaries
# # used as background grey reference lines in map in pdf
# nm_refBoundaries = here("_data","other_spatial","feature", "US_States.shp")  
# # background exclusion file (areas where you know presence points weren't collected [such as Tribal lands],
# # and thus where background points should be exluded)
# #   set up as geopackage, set as NULL if there are no exclusion areas
# #nm_bkgExclAreas <- c(here("_data","other_spatial","feature","az_tribal.gpkg"),"az_tribal")
# nm_bkgExclAreas <- NULL
# # background bias file (if presence points clearly show sampling bias, such as only near roads,
# # defining a distance-to raster here that defines the bias (e.g. distance to roads), will
# # subsample background points weighted by this bias)
# #  set to NULL if you don't want to define a bias file
# #nm_biasDistRas <- here("_data","other_spatial","raster","eucDist_AZ_roads.tif")
# nm_biasDistRas <- NULL
# 
# # project overview - this appears in the first paragraph of the metadata
# project_overview = "This model was developed for the Pennsylvania Natural Heritage Program."
# 
# # model comment in database
# model_comments = ""
# 
# # comment printed in PDF metadata
# metaData_comments = ""
# 
# 
# 
# 
# 
