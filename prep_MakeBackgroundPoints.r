library(sf)
library(here)
library(RSQLite)

studyArea_sf <- arc.open(studyArea)
studyArea_sf<- arc.select(studyArea_sf)
studyArea_sf <- arc.data2sf(studyArea_sf)

numpts <- 5000

# generate points
samps1 <- st_sample(studyArea_sf, size=numpts)
samps1 <- st_sf(fid = 1:length(samps1), geometry = samps1)
#samps <- st_join(samps1, stdyAreaHucs, join = st_intersects)[c("fid", "HUC10")]
#names(samps) <- c("fid", "huc10", "geometry")

st_crs(samps1) <- projPA



# get the predictor data ##############################################################################################
cat("Loading the predictor data...")
predictors_Current <- stack(list.files(here::here("_data","env_vars","Climate_Current"), pattern = 'asc$', full.names=TRUE ))

# extract values to point from the raster stack
bkgd_attributed <- raster::extract(predictors_Current, samps1, method="simple", sp=TRUE)
bkgd_attributed <- st_as_sf(points_attributed)