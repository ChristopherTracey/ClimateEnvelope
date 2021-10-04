library(stars)

# get the basemap data ################################################################################################
bndPA <- arc.open(here::here("modeling_data.gdb", "bnd_PAstate"))
bndPA <- arc.select(bndPA)
bndPA <- arc.data2sf(bndPA)

# get the predictor data ##############################################################################################
predictors_Current <- stack(list.files(here::here("_data","env_vars","Climate_Current"), pattern = 'asc$', full.names=TRUE ))
#predictors_Current <- read_stars(list.files(here::here("_data","env_vars","Climate_Current"), pattern='asc$', full.names=TRUE))
predictors_Future <- stack(list.files(here::here("_data","env_vars","Climate_Future2050rfp45"), pattern = 'asc$', full.names=TRUE ))
#predictors_Future <- read_stars(list.files(here::here("_data","env_vars","Climate_Future2050rfp45"), pattern='asc$', full.names=TRUE ))

# check to see if the names are the sames
setdiff(names(predictors_Current), names(predictors_Future))
setdiff(names(predictors_Future), names(predictors_Current))

# extract values to point from the raster stack
points_attributed <- extract(predictors_Current, species_sf, method="simple", sp=TRUE)
points_attributed <- st_as_sf(points_attributed)

library(vtable)
sumtable(points_attributed)

# look for duplicated values
tmp <-points_attributed
st_geometry(tmp) <- NULL
duplicated(tmp)










# set the modeling extent to be the extent of the rasters
ext = extent(predictors_Current)


