# install and/or load necessary packages and libraries
require(raster)
require(rgdal)
require(sf)
require(tidyverse)
require(fasterize)
require(here)
require(virtualspecies)
library(arcgisbinding)
library(RSQLite)
arc.check_product()

here::here()

####################
### Basic set up ###
####################

# set up projection parameter for use throughout script--albers
ascproj <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +") #projection of AdaptWest asc layers
#proj <- CRS("+proj=aea +lat_1=40 +lat_2=42 +lat_0=39 +lon_0=-78 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs") #projection for PA

# path to the shape to use for clipping--a shapefile in the same projection as rasters
studyArea <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/Refugia Modeling Boundary/Refugia Modeling Boundary.shp"


# get the basemap data ################################################################################################
studyArea <- arc.open(studyArea)
studyArea <- arc.select(studyArea)
studyArea <- arc.data2sf(studyArea)

#get extent of shape
clpextent<- st_bbox(studyArea)

#######################################
###### Clip a set of .tif rasters #####
#######################################

pathToTifs <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_NORM_8110_Bioclim_ASCII"
#"W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_NORM_8110_Bioclim_ASCII"
#W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_ENSEMBLE_rcp45_2050s_Bioclim_ASCII/
#W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_ENSEMBLE_rcp85_2050s_Bioclim_ASCII/

# the path to write out the clipped rasters to
pathToClipped <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_ENSEMBLE_rcp85_2050s_Bioclim_ASCII/NewVars_masked"
ifelse(!dir.exists(pathToClipped), dir.create(pathToClipped), FALSE)

# get a list of the grids, if asc already
tiflist <- list.files(path=pathToTifs, pattern=".tif$")

## already got some clipped? use the next few lines to check and 
## remove the ones already done
donetiflist <- list.files(path=pathToClipped, pattern=".tif$")
finaltifList <- tiflist[!tiflist %in% donetiflist]
tiflist <- finaltifList

# tack on the full paths and name them
gridlist <- as.list(paste(pathToTifs, tiflist, sep = "/"))
nm <- substr(tiflist, 1, nchar(tiflist) - 4)
names(gridlist) <- nm

# create list of full paths to write clipped rasters to
tifnm <- paste(nm, ".asc", sep="")
outfiles <- as.list(paste(pathToClipped, tifnm, sep= "/"))
names(outfiles) <- nm

## clip the rasters ----
for (i in 1:length(gridlist)){
  ras <- raster(gridlist[[i]], RAT=FALSE)
  fn <- paste(pathToClipped, "/", names(gridlist[i]), ".tif", sep="")
  a <- crop(ras, clpextent, filename=fn, format="GTiff", overwrite=TRUE)
}

predictors_Future <- stack(outfiles) #swapped this out manually to do the three sets of predictor variables

writeRaster(predictors_Future, names(predictors_Future), bylayer=TRUE, format="GTiff")


#stack raster layers and then examine for correlated layers
pathPredictorsCurrent <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_NORM_8110_Bioclim_ASCII/NewVars_masked"
pathPredictorsFuture4.5 <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_ENSEMBLE_rcp45_2050s_Bioclim_ASCII/NewVars_masked"
pathPredictorsFuture8.5 <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/EnvironmentalData_FullExtent/NA_ENSEMBLE_rcp85_2050s_Bioclim_ASCII/NewVars_masked"

# write file names into database here for the variables #


predictors_Current <- stack(list.files(pathPredictorsCurrent, pattern = 'tif$', full.names=TRUE )) #Reads in .tif files as a raster stack

###########################################################
# check for correlated variable groups and write to the database. This only needs to be done for the current climate variables
corr <- removeCollinearity(predictors_Current, plot=TRUE, multicollinearity.cutoff=0.8, select.variables=FALSE)
names(corr) <- seq(1:length(corr))
group <- unlist(mapply(function(x,y){ rep(y, length(x)) }, corr, names(corr)))
corr <- unlist(corr, use.names=TRUE)
corr <- as.data.frame(corr)
corr$group <- group


#plot(corr) #save image to folder w/ variables, for future reference

#this was overwriting in a weird order and I didn't feel like fixing automatically so I fixed manually. Come back and figure out why...
for(i in 1:nrow(corr)){
  db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
  dbSendQuery(db_cem, paste("UPDATE lkpEnvVar SET corrgroup = ", sQuote(corr$group[i])," WHERE rasname = ", sQuote(corr$corr[i]),";", sep="") )
  dbDisconnect(db_cem)
  #Sys.sleep(.1)
}



