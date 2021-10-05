#library(stars)
library(caret)

# get the basemap data ################################################################################################
projectArea <- arc.open(project_area)
projectArea <- arc.select(projectArea)
projectArea <- arc.data2sf(projectArea)

# get the predictor data ##############################################################################################
predictors_Current <- stack(list.files(here::here("_data","env_vars","Climate_Current"), pattern = 'asc$', full.names=TRUE ))
predictors_Future <- stack(list.files(here::here("_data","env_vars","Climate_Future2050rfp45"), pattern = 'asc$', full.names=TRUE ))
#predictors_Current <- read_stars(list.files(here::here("_data","env_vars","Climate_Current"), pattern='asc$', full.names=TRUE))
#predictors_Future <- read_stars(list.files(here::here("_data","env_vars","Climate_Future2050rfp45"), pattern='asc$', full.names=TRUE ))

# check to see if the names are the same
setdiff(names(predictors_Current), names(predictors_Future))
setdiff(names(predictors_Future), names(predictors_Current))

# extract values to point from the raster stack
points_attributed <- extract(predictors_Current, species_sf, method="simple", sp=TRUE)
points_attributed <- st_as_sf(points_attributed)


# look for near zero variation for each attributed variable ##############################################################
ptsAtt_df <-points_attributed
st_geometry(ptsAtt_df) <- NULL
ck_nrZeroVar <- nearZeroVar(ptsAtt_df, saveMetrics=TRUE)
if(any(ck_nrZeroVar$nzv==TRUE)){
  cat("Print at least one variable has near zero variation.")
} else {
  cat("There is variation is all the variables.")
}

# look for duplicated values ##############################################################
if(any(isTRUE(duplicated(tmp)))){
  cat("Print at least one variable has duplicated environmental variables.")
} else {
  cat("There are no fully duplicated variables.")
}

rm(tmp)








# set the modeling extent to be the extent of the rasters
ext = extent(predictors_Current)


