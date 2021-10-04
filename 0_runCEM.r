library(here)
library(dismo)
library(arcgisbinding)
arc.check_product()
library(ggplot2)

# species code (from lkpSpecies in modelling database. This will be the new folder name containing inputs/outputs)
sp_code <- "lupipere" # Lupinus perennis

# Modeling database
nm_db_file <- here("_data", "databases", "CEMdata.sqlite")




options(useFancyQuotes=FALSE) # needed to make sure SQL queries work as well as they could







#####


# loc_scripts is your repository. Make sure your git repository is set to correct branch
loc_scripts <- here()
# The main modelling folder for inputs/outputs. All sub-folders are created during the model run (when starting with step 1)
loc_model <- here("_data", "species")

# locations file (presence reaches). Provide full path; File is copied to modeling folder and timestamped.
nm_presFile <- here("_data", "occurrence", paste0(model_species, ".shp"))
#nm_presFile <- here("_data", "occurrence", paste0(sub("-","_",model_species), ".gpkg"))
# env vars location [Terrestrial-only variable]
loc_envVars = here("_data","env_vars","raster", "ras") # "D:/R_tmp/taenmont_20210304_134229"
# Name of background/envvars sqlite geodatabase, and base table name (2 length vector)
nm_bkgPts <- c(here("_data","env_vars","tabular", "background_CONUS.sqlite"), "background_pts")
# HUC spatial data set (shapefile) that is subsetted and used to define modeling area//range
nm_HUC_file <- here("_data","other_spatial","feature","HUC10.shp")
# map reference boundaries
# used as background grey reference lines in map in pdf
nm_refBoundaries = here("_data","other_spatial","feature", "US_States.shp")  
# background exclusion file (areas where you know presence points weren't collected [such as Tribal lands],
# and thus where background points should be exluded)
#   set up as geopackage, set as NULL if there are no exclusion areas
#nm_bkgExclAreas <- c(here("_data","other_spatial","feature","az_tribal.gpkg"),"az_tribal")
nm_bkgExclAreas <- NULL
# background bias file (if presence points clearly show sampling bias, such as only near roads,
# defining a distance-to raster here that defines the bias (e.g. distance to roads), will
# subsample background points weighted by this bias)
#  set to NULL if you don't want to define a bias file
#nm_biasDistRas <- here("_data","other_spatial","raster","eucDist_AZ_roads.tif")
nm_biasDistRas <- NULL

# project overview - this appears in the first paragraph of the metadata
project_overview = "This model was developed for the Pennsylvania Natural Heritage Program."

# model comment in database
model_comments = ""

# comment printed in PDF metadata
metaData_comments = ""

# your name
modeller = "Christopher Tracey"



source("1_PrepSpeciesData.r")
source("1_PrepSpeciesData.r")
source("1_PrepSpeciesData.r")