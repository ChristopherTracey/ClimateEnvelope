### Get the list of data sources

#source files

GBIF_data <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/GBIF_data_clip.shp"
iNat_data <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/iNat_data_clip.shp"
NatureServe_data <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/NSdata_cen.shp"
NY_plots1 <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/NY_plotdata_spplocations_XYTableToPoint.shp"
NY_plots2 <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/NY_Werier_data_XYTableToPoint.shp"
WV_Marshallia <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/WV_Marshallia_centroid.shp"
WV_plots <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/WV_Species_in_Plots_XYTableToPoint.shp"

GBIF_data <- arc.open(GBIF_data)
GBIF_data <- arc.select(GBIF_data)
GBIF_data <- arc.data2sf(GBIF_data)

NatureServe_data <- arc.open(NatureServe_data)
NatureServe_data <- arc.select(NatureServe_data)
NatureServe_data <- arc.data2sf(NatureServe_data)
