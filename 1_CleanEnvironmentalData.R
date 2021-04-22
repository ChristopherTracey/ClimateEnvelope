# install and/or load necessary packages and libraries

if (!requireNamespace("raster", quietly = TRUE)) install.packages("raster")
require(raster)

if (!requireNamespace("rgdal", quietly = TRUE)) install.packages("rgdal")
require(rgdal)

if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
require(sf)

if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
require(tidyverse)

if (!requireNamespace("fasterize", quietly = TRUE)) install.packages("fasterize")
require(fasterize)

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
require(here)

if (!requireNamespace("virtualspecies", quietly = TRUE)) install.packages("virtualspecies")
require(virtualspecies)

here::i_am("Scripting/1_CleanEnvironmentalData.R")

####################
### Basic set up ###
####################


# set up projection parameter for use throughout script--albers
ascproj <- CRS("+proj=lcc +lon_0=-95 +lat_1=49 +lat_2=77 +lat_0=0") #projection of AdaptWest asc layers
proj <- CRS("+proj=aea +lat_1=40 +lat_2=42 +lat_0=39 +lon_0=-78 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs") #projection for PA

# path to the shape to use for clipping--a shapefile in the same projection as rasters
pathToClipShape <- "C:/Users/AJohnson/Documents/ArcGIS/Projects/ClimateRefugia/EnvironmentalData"
clipShapeName <- "States_BoundingBox"

clpShp <- readOGR(pathToClipShape,clipShapeName)

#Reproject the clip-shape to match the asc layers
reproj_clpshp <- spTransform(clpShp, CRS=ascproj)

#get extent of shape
clpextent<- extent(reproj_clpshp)


#######################################
###### Clip a set of .asc rasters #####
#######################################

pathToASCs <- "C:/Users/AJohnson/Documents/ArcGIS/Projects/ClimateRefugia/EnvironmentalData_FullExtent/NA_ENSEMBLE_rcp45_2050s_Bioclim_ASCII/"

# the path to write out the clipped rasters to
pathToClipped <- "C:/Users/AJohnson/Documents/ArcGIS/Projects/ClimateRefugia/EnvironmentalData_FullExtent/NA_ENSEMBLE_rcp45_2050s_Bioclim_ASCII/NewVars_masked/"

# get a list of the grids, if asc already
asclist <- list.files(path = pathToASCs, pattern = ".asc$")

## already got some clipped? use the next few lines to check and 
## remove the ones already done
doneasclist <- list.files(path = pathToClipped, pattern = ".asc$")
finalascList <- asclist[!asclist %in% doneasclist]
asclist <- finalascList

# tack on the full paths and name them
gridlist<-as.list(paste(pathToASCs, asclist,sep = "/"))
nm <- substr(asclist,1,nchar(asclist) - 4)
names(gridlist)<-nm

# create list of full paths to write clipped rasters to
ascnm <- paste(nm, ".asc", sep="")
outfiles <- as.list(paste(pathToClipped, ascnm, sep= "/"))
names(outfiles)<-nm

## clip the .asc rasters ----
for (i in 1:length(gridlist)){
  ras <- raster(gridlist[[i]])
  fn <- paste(pathToClipped, "/", names(gridlist[i]), ".asc", sep="")
  a <- crop(ras,reproj_clpshp, filename = fn, format = "ascii", overwrite = TRUE)
}

####################################################
##### Convert a set of .tif files to .asc files ####
####################################################

#path to the folder where the list of .tif files lives
pathToTIFs <- "C:/Users/AJohnson/Documents/ArcGIS/Projects/ClimateRefugia/EnvironmentalData_FullExtent/LandscapeCondition"

# get a list of the .tif grids
tiflist <- list.files(path = pathToTIFs, pattern = ".tif$")

#create list of full paths to convert to .asc files
gridlist<-as.list(paste(pathToTIFs, tiflist,sep = "/"))
nm <- substr(tiflist,1,nchar(tiflist) - 4)
nm <- paste(nm, ".asc", sep= "")
names(gridlist)<-nm

#turn into a raster and write out as an .asc file
ras <- list()
for (i in 1:length(gridlist)){
  ras[i] <- raster(gridlist[[i]])
  writeRaster(ras[[i]], filename=paste(pathToTIFs, nm[[i]]), format="ascii", overwrite=TRUE)
}

#######################################################
#### Create raster stack from clipped raster layers ###
#######################################################

# tack on the full paths and name them
cliplist<-as.list(paste(pathToClipped, doneasclist,sep = "/"))
nm <- substr(doneasclist,1,nchar(doneasclist) - 4)
names(cliplist)<-nm


#stack raster layers and then examine for correlated layers
envtStack <- stack(cliplist) #Reads in .asc files as a raster stack

corr <- removeCollinearity(envtStack, plot = TRUE, multicollinearity.cutoff = 0.8)
plot(corr) #save image to folder w/ variables, for future reference


# Automatic selection of variables not intercorrelated (likely doesn't choose the "right" reps)
uncorr <- removeCollinearity(envtStack, plot = TRUE, multicollinearity.cutoff = 0.8, select.variables = TRUE)
length(uncorr)

#custom list of length uncorr of variable names
UncorrVars <- c("MAR","RH","MAP","MSP","CMD","TD","DD5","PAS")
UncorrVars <- paste(UncorrVars,".asc", sep="")

#path to uncorrelated variable folder
pathToUnCorr <- "C:/Users/AJohnson/Documents/ArcGIS/Projects/ClimateRefugia/EnvironmentalData/NA_NORM_8110_Bioclim_ASCII/NewVars_masked/UnCorr"


uncorrlist<-as.list(paste(pathToUnCorr, UncorrVars,sep = "/"))
nm <- substr(UncorrVars,1,nchar(UncorrVars) - 4)
names(uncorrlist)<-nm

fn <- 
for (i in 1:length(uncorrlist)){
  a <- writeRaster(uncorrlist, filename = uncorrlist, format = "ascii", overwrite = TRUE)
}


