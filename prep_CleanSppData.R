### Get the list of data sources

#source files

#Species list
SpList <- read.csv(file="W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/SpeciesList_DataTracking.csv")

SpList_modelready <- SpList[(SpList$Model.Ready=="Y"),] #just the species ready to model (i.e. we have enough point data ready to go)
ELCODEs <- SpList %>% dplyr::select(ELCODE, SNAME) #list of all the ELCODES in this modeling project w/ SNAME as key column to join to the species point data files

#Species point files
Biotics_pointreps <-"W:/Heritage/Heritage_Data/Biotics_datasets.gdb/eo_ptreps"
GBIF_data <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/GBIF_data_clip.shp"
iNat_data <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/iNat_data_clip.shp"
NatureServe_data <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/NSdata_cen.shp"
NY_plots1 <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/NY_plotdata_spplocations_XYTableToPoint.shp"
NY_plots2 <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/NY_Werier_data_XYTableToPoint.shp"
WV_Marshallia <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/WV_Marshallia_centroid.shp"
WV_plots <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/WV_Species_in_Plots_XYTableToPoint.shp"
PA_plots <- "W:/Heritage/Heritage_Projects/1280_CC_Refugia/SppSpatialData/PA_Species_in_plots.shp"

#reformat each sf so they can be joined into a single dataset for modeling

# Columns for the final dataset
SpNames <- c("ELCODE","SNAME","YEAR","NOTES","DATASOURCE","USEDATA", "geom")

CRS_Code <- 5070 #NAD83 Conus Albers

###### Biotics ##########
pointreps <- arc.open(Biotics_pointreps)

SQLquery_pointreps <- paste("'LASTOBS_YR > 1991'")
selected_pointreps <- arc.select(pointreps, c('ELCODE', 'SNAME', 'ID_CONFIRM','LASTOBS_YR', 'EORANK', 'EO_DATA'), where_clause="LASTOBS_YR > 1990 AND EORANK NOT IN ('X','F') AND ID_CONFIRM = 'Y - Yes'")
selected_pointreps <- arc.data2sf(selected_pointreps)
selected_pointreps <- filter(selected_pointreps, ELCODE %in% SpList$ELCODE)
selected_pointreps$datasource <- "PA Biotics"
selected_pointreps$usedata <- "Y"

selected_pointreps <- selected_pointreps %>% dplyr::select(ELCODE, SNAME, LASTOBS_YR, EO_DATA, datasource, usedata)
names(selected_pointreps) <- SpNames

str(selected_pointreps)
selected_pointreps <- st_transform(selected_pointreps, crs=CRS_Code) #the CRS code for Biotics

###### GBIF ##############
GBIF_data <- arc.open(GBIF_data)
GBIF_data <- arc.select(GBIF_data)
GBIF_data <- arc.data2sf(GBIF_data)

GBIF_data <- right_join(GBIF_data, ELCODEs)
GBIF_data <- GBIF_data %>% dplyr::select(ELCODE, SNAME,LastObs, Notes, DataSource, usedata)
names(GBIF_data) <- SpNames

GBIF_data <- st_transform(GBIF_data, crs=CRS_Code)

###### iNat Data ############
iNat_data <- arc.open(iNat_data)
iNat_data <- arc.select(iNat_data)
iNat_data <- arc.data2sf(iNat_data)
iNat_data$Notes <- NA

iNat_data <- right_join(iNat_data, ELCODEs)
iNat_data <- iNat_data %>% dplyr::select(ELCODE, SNAME,LastObs, Notes, DataSource, usedata)
names(iNat_data) <- SpNames

iNat_data <- st_transform(iNat_data, crs=CRS_Code)


###### 2 NatureServe #########
NatureServe_data <- arc.open(NatureServe_data)
NatureServe_data <- arc.select(NatureServe_data)
NatureServe_data <- arc.data2sf(NatureServe_data)

NatureServe_data$Notes <- NA
NatureServe_data$usedata <- "Y"

NatureServe_data <- right_join(NatureServe_data, ELCODEs)
NatureServe_data <- NatureServe_data %>% dplyr::select(ELCODE, SNAME, LASTOBS_YR, Notes, source, usedata)
names(NatureServe_data) <- SpNames

NatureServe_data <- NatureServe_data[(NatureServe_data$YEAR >1990),] #turns out there are a lot of old records in here
NatureServe_data <- st_transform(NatureServe_data, crs=CRS_Code)

###### NY plots 1 #########
NY_plots1 <- arc.open(NY_plots1)
NY_plots1 <- arc.select(NY_plots1)
NY_plots1 <- arc.data2sf(NY_plots1)
NY_plots1$usedata <- "Y"

NY_plots1 <- NY_plots1 %>% dplyr::select(ELCODE, SNAME, surveydate, comments, datasource, usedata)
names(NY_plots1) <- SpNames

NY_plots1$YEAR <- year(NY_plots1$YEAR) #extract just the year

NY_plots1$DATASOURCE <- "NY Plots data"

NY_plots1 <- st_transform(NY_plots1, crs=CRS_Code)

###### NY plots 2 #########
NY_plots2 <- arc.open(NY_plots2)
NY_plots2 <- arc.select(NY_plots2)
NY_plots2 <- arc.data2sf(NY_plots2)
NY_plots2$usedata <- "Y"
NY_plots2$Notes <- NA

NY_plots2 <- NY_plots2 %>% dplyr::select(ELCODE, SNAME, surveydate, Notes, datasource, usedata)
names(NY_plots2) <- SpNames

NY_plots2$YEAR <- NA

NY_plots2 <- st_transform(NY_plots2, crs=CRS_Code)

####### Marshallia ########
WV_Marshallia <- arc.open(WV_Marshallia)
WV_Marshallia <- arc.select(WV_Marshallia)
WV_Marshallia <- arc.data2sf(WV_Marshallia)

WV_Marshallia$datasource <- "WV Marshallia EOs"
WV_Marshallia$usedata <- "Y"

WV_Marshallia <- WV_Marshallia %>% dplyr::select(ELCODE, SNAME,LAST_OBS_D, BASIC_EO_R, datasource, usedata)
names(WV_Marshallia) <- SpNames

WV_Marshallia$YEAR <- ymd(WV_Marshallia$YEAR)
WV_Marshallia$YEAR <- year(WV_Marshallia$YEAR)

WV_Marshallia <- st_transform(WV_Marshallia, crs=CRS_Code)

####### WV_plots #########
WV_plots <- arc.open(WV_plots)
WV_plots <- arc.select(WV_plots)
WV_plots <- arc.data2sf(WV_plots)

WV_plots$datasource <- "WV plot data"
WV_plots$usedata <- "Y"
WV_plots$Notes <- NA

WV_plots <- WV_plots %>% dplyr::select(ELCODE, SNAME, surveydate, Notes, datasource, usedata)
names(WV_plots) <- SpNames

WV_plots$YEAR <- NA

WV_plots <- st_transform(WV_plots, crs=CRS_Code)

####### PA_plots #########
PA_plots <- arc.open(PA_plots)
PA_plots <- arc.select(PA_plots)
PA_plots <- arc.data2sf(PA_plots)

PA_plots$usedata <- "Y"
PA_plots$Notes <- NA
PA_plots$surveydate <- NA
PA_plots$datasource <- "PA plot data"

PA_plots <- right_join(PA_plots, ELCODEs)
PA_plots <- PA_plots %>% dplyr::select(ELCODE, SNAME, surveydate, Notes, datasource, usedata)
names(PA_plots) <- SpNames

PA_plots <- st_transform(PA_plots, crs=CRS_Code)

#############################
#### Bind em all together ###
#############################

Sp_DataList <- list(selected_pointreps, GBIF_data, iNat_data, NatureServe_data, NY_plots1, NY_plots2, WV_Marshallia, WV_plots, PA_plots)

Sp_Points <- bind_rows(Sp_DataList)
Sp_Points %>% drop_na(geom)

st_write(Sp_Points, dsn = "Sp_Points.shp", layer = "Sp_Points.shp", driver = "ESRI Shapefile", append=FALSE)

ggplot() + geom_sf(data=Sp_Points, aes(color=DATASOURCE))
