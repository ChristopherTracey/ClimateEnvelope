



# get the basemap data ################################################################################################
bndPA <- arc.open(here::here("modeling_data.gdb", "bnd_PAstate"))
bndPA <- arc.select(bndPA)
bndPA <- arc.data2sf(bndPA)

# get the predictor data ##############################################################################################
predictors_Current <- stack(list.files(here::here("_data","env_vars","Climate_Current"), pattern = 'asc$', full.names=TRUE ))
predictors_Future <- stack(list.files(here::here("_data","env_vars","Climate_Future2050rfp45"), pattern = 'asc$', full.names=TRUE ))
# check to see if the names are the sames
setdiff(names(predictors_Current), names(predictors_Future))
setdiff(names(predictors_Future), names(predictors_Current))

# set the modeling extent to be the extent of the rasters
ext = extent(predictors_Current)
