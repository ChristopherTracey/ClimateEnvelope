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
studyArea <- here::here("_data","other_spatial","modeling_data.gdb", "bound_pro")

# get the basemap data ################################################################################################
studyArea <- arc.open(studyArea)
studyArea <- arc.select(studyArea)
studyArea <- arc.data2sf(studyArea)

#get extent of shape
clpextent<- st_bbox(studyArea)

#######################################
###### Clip a set of .tif rasters #####
#######################################

pathToTifs <- "E:/Refugia/climateData/ensemble_ssp245_2011_bioclim"

# the path to write out the clipped rasters to
pathToClipped <-  here::here("_data","env_vars",basename(pathToTifs))
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


# rename files
srcString <- str_replace(basename(pathToTifs), "bioclim", "")

fromfiles <- list.files(pathToClipped, pattern=srcString)
tofiles <- str_replace(fromfiles, srcString, "")

file.rename(file.path(pathToClipped, fromfiles), file.path(pathToClipped, tofiles))

# tack on the full paths and name them
cliplist <- as.list(paste(pathToClipped, tofiles, sep="/")) 
nm <- substr(tofiles, 1, nchar(tofiles) - 4)
names(cliplist) <- nm

#stack raster layers and then examine for correlated layers
predictors_Current <- stack(cliplist) #Reads in .asc files as a raster stack

###########################################################
# check for correlated variable groups and write to the database
corr <- removeCollinearity(predictors_Current, plot=TRUE, multicollinearity.cutoff=0.8, select.variables=FALSE)
names(corr) <- seq(1:length(corr))
corr <- unlist(corr)
corr <- as.data.frame(corr)
corr$group <- rownames(corr)
corr$group <- substr(corr$group, 1, 1)
corr <- corr[corr$group %in% corr$group[duplicated(corr$group)],]

#plot(corr) #save image to folder w/ variables, for future reference

for(i in 1:nrow(corr)){
  db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
  dbSendQuery(db_cem, paste("UPDATE lkpEnvVar SET CorrGroup = ", sQuote(corr$group[i])," WHERE rasCode = ", sQuote(corr$corr[i]),";", sep="") )
  dbDisconnect(db_cem)
  #Sys.sleep(.1)
}
