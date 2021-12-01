# install and/or load necessary packages and libraries
require(sf)
require(tidyverse)
require(here)
library(arcgisbinding)
library(RSQLite)
arc.check_product()

here::here()

# get the list of species we're working on and build sqlite table
path_splist <- "P:/Conservation Programs/Natural Heritage Program/Projects/Active projects/1280 WRCP CC Refugia/Project Materials and Data/Species Lists/SpeciesList_CC_Refugia_Fall2021.csv"

splist <- read.csv(path_splist, stringsAsFactors=FALSE)
splist$cutecode <- paste(tolower(substr(word(splist$SNAME, 1), 1, 4)), substr(word(splist$SNAME, 2), 1, 4), sep="")  
length(unique(splist$cutecode))==nrow(splist)

splist_dup <- splist[duplicated(splist$cutecode),]   #subset(splist, duplicated(splist$cutecode))
splist <- splist[!duplicated(splist$cutecode),]
splist_dup$cutecode <- paste(substr(splist_dup$cutecode, 1, 7), 1, sep="")
splist <- rbind(splist, splist_dup)
splist <- splist[order(splist$SNAME),]
rm(splist_dup)

db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
# delete existing threats and recs for this site if they exist
dbExecute(db_cem, paste("DELETE FROM lkpSpecies"))
# add in the new data
dbAppendTable(db_cem, "lkpSpecies", splist)
dbDisconnect(db_cem)

########################
# get data
sp_PABiotics <- "W:/Heritage/Heritage_Data/Biotics_datasets.gdb" 
sp_NSbld <- "PA_CC_Refugia_EO_Polys_112021.shp"
sp_iNat <- "iNat_data_clip.shp"
sp_GBIF <- "GBIF_data_clip.shp"  

# get PA Biotics Data
ptreps <- arc.open(paste(sp_PABiotics, "eo_ptreps", sep="/"))  
ptreps <- arc.select(ptreps, c("SNAME","EO_ID","EST_RA","PREC_BCD","LASTOBS_YR"), where_clause=paste("SNAME IN (", paste(sQuote(splist$SNAME), collapse=", "), ")"))
ptreps <- arc.data2sf(ptreps)
ptreps$SUBNATION <- "PA"
ptreps$source <- "PA Biotics"

setdiff(ptreps$SNAME, splist$SNAME)
setdiff(splist$SNAME, ptreps$SNAME)

#ptreps <- st_transform(ptreps, ascproj) # reproject data

# get NS data
NSdata <- arc.open(here::here("_data","occurrence",sp_NSbld))
NSdata <- arc.select(NSdata, fields=c("SUBNATION","EO_ID","SNAME","LOBS_Y", "PA_SppList", "EST_REP_AC", "PREC_BCD"))
NSdata <- arc.data2sf(NSdata)
names(NSdata)[names(NSdata)=="LOBS_Y"] <- "LASTOBS_YR"
names(NSdata)[names(NSdata)=="EST_REP_AC"] <- "EST_RA"
NSdata$source <- "NatureServe BLD"

setdiff(names(NSdata), names(ptreps))
setdiff(names(ptreps), names(NSdata))

# find related infrataxa
a <- unique(st_drop_geometry(NSdata[which(NSdata$PA_SppList=="Y-Rel Infrataxa"),"SNAME"]))
a <- a$SNAME
a <- sort(a)

NSdata$SNAME[NSdata$SNAME=="Andromeda polifolia var. glaucophylla"] <- "Andromeda polifolia"            
NSdata$SNAME[NSdata$SNAME=="Andropogon glomeratus var. glomeratus"] <- "Andropogon glomeratus"            
NSdata$SNAME[NSdata$SNAME=="Aristida purpurascens var. purpurascens"] <- "Aristida purpurascens"         
NSdata$SNAME[NSdata$SNAME=="Aristida virgata"] <- "Aristida purpurascens"                                 
NSdata$SNAME[NSdata$SNAME=="Baptisia australis var. australis"] <- "Baptisia australis"                
NSdata$SNAME[NSdata$SNAME=="Bartonia paniculata ssp. paniculata"] <- "Bartonia paniculata"             
NSdata$SNAME[NSdata$SNAME=="Boltonia asteroides var. asteroides"] <- "Boltonia asteroides"              
NSdata$SNAME[NSdata$SNAME=="Boltonia asteroides var. glastifolia"] <- "Boltonia asteroides"             
NSdata$SNAME[NSdata$SNAME=="Bouteloua curtipendula var. curtipendula"] <- "Bouteloua curtipendula"        
NSdata$SNAME[NSdata$SNAME=="Carex lasiocarpa var. americana"] <- "Carex lasiocarpa"                  
NSdata$SNAME[NSdata$SNAME=="Carex oligosperma var. oligosperma"] <- "Carex oligosperma"
NSdata$SNAME[NSdata$SNAME=="Carex tetanica var. canbyi"] <- "Carex tetanica"
NSdata$SNAME[NSdata$SNAME=="Dichanthelium oligosanthes var. oligosanthes"] <- "Dichanthelium oligosanthes"
NSdata$SNAME[NSdata$SNAME=="Dichanthelium oligosanthes var. scribnerianum"] <- "Dichanthelium oligosanthes"
NSdata$SNAME[NSdata$SNAME=="Hypericum densiflorum var. interior"] <- "Hypericum densiflorum"
NSdata$SNAME[NSdata$SNAME=="Liatris scariosa var. nieuwlandii"] <- "Liatris scariosa"
NSdata$SNAME[NSdata$SNAME=="Liatris scariosa var. novae-angliae"] <- "Liatris scariosa"
NSdata$SNAME[NSdata$SNAME=="Linnaea borealis ssp. americana"] <- "Linnaea borealis"                 
NSdata$SNAME[NSdata$SNAME=="Lupinus perennis ssp. perennis"] <- "Lupinus perennis"
NSdata$SNAME[NSdata$SNAME=="Platanthera blephariglottis var. blephariglottis"] <- "Platanthera blephariglottis"
NSdata$SNAME[NSdata$SNAME=="Polygala cruciata var. aquilonia"] <- "Polygala cruciata"                
NSdata$SNAME[NSdata$SNAME=="Prunus alleghaniensis var. alleghaniensis"] <- "Prunus alleghaniensis"
NSdata$SNAME[NSdata$SNAME=="Ranunculus pusillus var. pusillus"] <- "Ranunculus pusillus"
NSdata$SNAME[NSdata$SNAME=="Rubus pubescens var. pubescens"] <- "Rubus pubescens"                  
NSdata$SNAME[NSdata$SNAME=="Rudbeckia fulgida var. fulgida"] <- "Rudbeckia fulgida"
NSdata$SNAME[NSdata$SNAME=="Scheuchzeria palustris ssp. americana"] <- "Scheuchzeria palustris"
NSdata$SNAME[NSdata$SNAME=="Solidago rigida var. glabrata"] <- "Solidago rigida"                   
NSdata$SNAME[NSdata$SNAME=="Solidago rigida var. rigida"] <- "Solidago rigida"
NSdata$SNAME[NSdata$SNAME=="Solidago uliginosa var. uliginosa"] <- "Solidago uliginosa"
NSdata$SNAME[NSdata$SNAME=="Spiranthes casei var. casei"] <- "Spiranthes casei"                     
NSdata$SNAME[NSdata$SNAME=="Stellaria borealis ssp. borealis"] <- "Stellaria borealis"  

# more synomomy
NSdata$SNAME[NSdata$SNAME=="Eupatorium coelestinum"] <- "Conoclinium coelestinum" 
NSdata$SNAME[NSdata$SNAME=="Sibbaldia tridentata"] <- "Potentilla tridentata" 
NSdata$SNAME[NSdata$SNAME=="Sibbaldiopsis tridentata"] <- "Potentilla tridentata" 
NSdata$SNAME[NSdata$SNAME=="Prunus susquehanae"] <- "Prunus pumila var. susquehanae" 
NSdata$SNAME[NSdata$SNAME=="Juncus balticus var. littoralis"] <- "Juncus arcticus var. littoralis" 
NSdata$SNAME[NSdata$SNAME=="Ripariosida hermaphrodita"] <- "Sida hermaphrodita" 
NSdata$SNAME[NSdata$SNAME=="Viburnum opulus var. americanum"] <- "Viburnum trilobum" 
NSdata$SNAME[NSdata$SNAME=="Hypericum ascyron ssp. pyramidatum"] <- "Hypericum pyramidatum" 
NSdata$SNAME[NSdata$SNAME=="Talinum teretifolium"] <- "Hypericum pyramidatum" 
NSdata$SNAME[NSdata$SNAME=="Patis racemosa"] <- "Piptatherum racemosum" 
NSdata$SNAME[NSdata$SNAME=="Arnoglossum muehlenbergii"] <- "Arnoglossum reniforme" 
NSdata$SNAME[NSdata$SNAME=="Borodinia missouriensis"] <- "Arabis missouriensis" 
NSdata$SNAME[NSdata$SNAME=="Phlox latifolia"] <- "Phlox ovata" 
NSdata$SNAME[NSdata$SNAME=="Prunus pumila var. cuneata"] <- "Prunus pumila var. susquehanae" 
NSdata$SNAME[NSdata$SNAME=="Crocanthemum propinquum"] <- "Helianthemum propinquum" 

setdiff(NSdata$SNAME, splist$SNAME)
setdiff(splist$SNAME, NSdata$SNAME)

NSdata$PA_SppList <- NULL

# project the data into a common projections

NSdata = st_cast(NSdata, "POLYGON") # i don't know why casting this to polygon fixes the crqazy projection problem, but it does.
st_write(NSdata, "NSdata.shp", delete_layer=TRUE)
NSdata <- st_transform(NSdata, st_crs(ptreps)) # reproject data to a common system
st_write(NSdata, "NSdata.shp", delete_layer=TRUE)


#NSdata_centroid <- st_centroid(NSdata)
NSdata_point <- st_point_on_surface(NSdata)
st_write(NSdata_point, "NSdata_cen.shp", delete_layer=TRUE)


#####################################
# get iNat data together

inat <- arc.open(here::here("_data","occurrence",sp_iNat))  
inat <- arc.select(inat)
inat <- arc.data2sf(inat)
inat$SUBNATION <- ""
inat$EST_RA <- ""
inat$PREC_BCD <- ""
names(inat)[names(inat)=="LastObs"] <- "LASTOBS_YR"
names(inat)[names(inat)=="DataSource"] <- "source"
names(inat)[names(inat)=="DataID"] <- "EO_ID"
inat$EO_ID <- str_replace(inat$EO_ID, "https://www.inaturalist.org/observations/", "")
inat$FID <- NULL

setdiff(inat$SNAME, splist$SNAME)
setdiff(splist$SNAME, inat$SNAME)

inat <- st_transform(inat, st_crs(ptreps))

inatDataSummary <- as.data.frame(table(inat$SNAME))
names(inatDataSummary)[names(inatDataSummary)=="Var1"] <- "SNAME"
names(inatDataSummary)[names(inatDataSummary)=="Freq"] <- "Record_inat"

#####################################
# get GBIF data together

gbif <- arc.open(here::here("_data","occurrence",sp_GBIF))  
gbif <- arc.select(gbif)
gbif <- arc.data2sf(gbif)
gbif$SUBNATION <- ""
gbif$EST_RA <- ""
gbif$PREC_BCD <- ""
names(gbif)[names(gbif)=="LastObs"] <- "LASTOBS_YR"
names(gbif)[names(gbif)=="DataSource"] <- "source"
names(gbif)[names(gbif)=="DataID"] <- "EO_ID"
gbif$FID <- NULL
gbif$Notes <- NULL
gbif$usedata <- NULL
gbif$coordinate <- NULL

gbif <- st_transform(gbif, st_crs(ptreps))

setdiff(gbif$SNAME, splist$SNAME)
setdiff(splist$SNAME, gbif$SNAME)

# create summary of the gbif data
gbifDataSummary <- as.data.frame(table(gbif$SNAME))
names(gbifDataSummary)[names(gbifDataSummary)=="Var1"] <- "SNAME"
names(gbifDataSummary)[names(gbifDataSummary)=="Freq"] <- "Record_gbif"

####################################
# join the EO layers together
spData <- rbind(ptreps, NSdata_point, gbif, inat)

spData <- merge(spData, splist, by="SNAME", all.x=TRUE)


#######################################
# do some data cleaning
tapply(spData$LASTOBS_YR, spData$SUBNATION, summary)# summary of last obs year by state

# create summary of the original data
spDataSummary <- as.data.frame.matrix(table(spData$SNAME, spData$source))
spDataSummary <- data.frame(SNAME=row.names(spDataSummary), spDataSummary)
rownames(spDataSummary) <- NULL
names(spDataSummary)[names(spDataSummary)=="NatureServe.BLD"] <- "NatureServeBLD"
names(spDataSummary)[names(spDataSummary)=="PA.Biotics"] <- "PABiotics"
# remove data beyond our cutoff year
cutoffYear <- 1991
spData <- spData[which(spData$LASTOBS_YR>=cutoffYear),]
# add the data lost by the cutoffyear to the summary table
tmp_spDataSummary <- as.data.frame.matrix(table(spData$SNAME, spData$source))
tmp_spDataSummary<- data.frame(SNAME=row.names(tmp_spDataSummary), tmp_spDataSummary)
rownames(tmp_spDataSummary) <- NULL
tmp_spDataSummary$GBIF <- NULL
tmp_spDataSummary$iNaturalist <- NULL
names(tmp_spDataSummary)[names(tmp_spDataSummary)=="NatureServe.BLD"] <- "NatureServeBLD_1991"
names(tmp_spDataSummary)[names(tmp_spDataSummary)=="PA.Biotics"] <- "PABiotics1991"
spDataSummary <- merge(spDataSummary, tmp_spDataSummary, by="SNAME", all.x=TRUE)
rm(tmp_spDataSummary)

spDataSummary$total1991 <- spDataSummary$GBIF + spDataSummary$iNaturalist + spDataSummary$NatureServeBLD_1991 + spDataSummary$PABiotics1991

# write a master copy of the data
st_write(spData, "SpeciesDataMaster.shp", delete_layer=TRUE)

#generate individual maps for training
looplist <- unique(spData$cutecode)
looplist <- looplist[!is.na(looplist)]

############
StateBnd <- arc.open(here::here("_data","other_spatial","RefugiaModelingBoundary.shp"))  
StateBnd <- arc.select(StateBnd)
StateBnd <- arc.data2sf(StateBnd)
StateBnd <- st_transform(StateBnd, st_crs(ptreps))

######################
# BONAP

bonap <- read.csv("W:/Heritage/Heritage_Projects/1495_PlantConservationPlan/BONAP/bonap_data.csv", stringsAsFactors=FALSE)

bonap1 <- pivot_longer(bonap, cols=10:1940)
bonap1$name <- str_replace(bonap1$name, "_", " ") 
bonap1$name <- str_replace(bonap1$name, " var ", " var. ")
bonap1$name <- str_replace(bonap1$name, " ssp ", " ssp. ") 

bonap2 <- bonap1[which(bonap1$name %in% splist$SNAME),]
bonap2 <- bonap2[which(bonap2$STATEFP %in% c(10,11,24,34,36,39,42,51,54)),]

counties <- arc.open(here::here("_data","other_spatial","studyArea_counties.shp"))  
counties <- arc.select(counties)
counties <- arc.data2sf(counties)
counties <- st_transform(counties, st_crs(ptreps))

bonap2$FIPS <- paste0(str_pad(bonap2$STATEFP,2,side="left", pad="0"), str_pad(bonap2$COUNTYFP,3,side="left", pad="0"))
bonap2 <- merge(bonap2, splist[c("SNAME","cutecode")], by.y="SNAME", by.x="name", all.x=TRUE)
bonap2 <- bonap2[which(bonap2$value=="PNR"|bonap2$value=="PR"|bonap2$value=="EXT"),]


#######################

spData$source <- ordered(spData$source, levels=c("PA Biotics","NatureServe BLD","GBIF","iNaturalist"))

for(i in 1:length(looplist)){ # 
  # subset the data frame for the species
  spData_sub <- spData[which(spData$cutecode==looplist[i]),]
  # get bonap data
  bonap_sub <- bonap2[which(bonap2$cutecode==looplist[i]),]
  bonap_sub <- merge(counties, bonap_sub, by.x="COUNTY_KEY", by.y="FIPS")
    
  # make the map
  p <- ggplot() +
    geom_sf(data=bonap_sub, aes(), colour=NA, fill="gray") +
    geom_sf(data=StateBnd, aes(), colour="black", fill=NA) +
    geom_sf(data=spData_sub, mapping=aes(color=source), alpha=0.9) + # +
    scale_color_manual(
      breaks=c("PA Biotics","NatureServe BLD","GBIF","iNaturalist"), 
      values=c("PA Biotics"="#2471a3", "NatureServe BLD"="#16a085", "GBIF"="#f39c12", "iNaturalist"="#b03a2e"),
      labels=c("PA Biotics","NatureServe BLD","GBIF","iNaturalist"), drop=FALSE, na.value="white") + #
    labs(
      title=unique(spData_sub$SNAME)
    ) +
    theme_void() +
    theme(legend.position=c(0.1,0.85)) +
    theme(legend.title=element_blank()) +
    theme(legend.text=element_text(size=10)) +
    theme(axis.text=element_blank(), axis.title=element_text(size=15)) 
   
    ggsave(filename=paste(here::here("_data","other_spatial","speciesmaps"),"/","spmap_",looplist[i],".png",sep=""), plot=p, width=14, height=10, units="cm", scale=1, dpi=300)
}




