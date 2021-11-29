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


sp_PABiotics <- "" 
sp_NSbld <- "PA_CC_Refugia_EO_Polys_112021.shp"
sp_iNat <- NULL
sp_GBIF <- NULL  


# get the study area data ################################################################################################
studyArea <- arc.open(studyArea)
studyArea <- arc.select(studyArea)
studyArea <- arc.data2sf(studyArea)




