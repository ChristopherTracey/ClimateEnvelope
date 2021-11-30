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
sp_iNat <- NULL
sp_GBIF <- NULL  

# get PA Biotics Data
ptreps <- arc.open(paste(sp_PABiotics, "eo_ptreps", sep="/"))  
ptreps <- arc.select(ptreps, c("SNAME","SCOMNAME","EO_ID","EST_RA","PREC_BCD","LASTOBS_YR"), where_clause=paste("SNAME IN (", paste(sQuote(splist$SNAME), collapse=", "), ")"))
ptreps <- arc.data2sf(ptreps)
ptreps$SUBNATION <- "PA"
ptreps$source <- "PA Biotics"

setdiff(ptreps$SNAME, splist$SNAME)
setdiff(splist$SNAME, ptreps$SNAME)

ptreps <- st_transform(ptreps, ascproj) # reproject data

# get NS data
NSdata <- arc.open(here::here("_data","occurrence",sp_NSbld))
NSdata <- arc.select(NSdata, fields=c("SUBNATION","EO_ID","SNAME","SCOMNAME","LOBS_Y", "PA_SppList", "EST_REP_AC", "PREC_BCD"))
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

NSdataSinglepart = st_cast(NSdata, "POLYGON") # i don't know why casting this to polygon fixes the crqazy projection problem, but it does.

NSdata <- st_transform(NSdata, st_crs(ptreps)) # reproject data to a common system
st_write(NSdata, "NSdata.shp")
NSdata_centroid <- st_centroid(NSdata)
st_write(NSdata_centroid, "NSdata_cen.shp")

####################################
# join the layers together
spData <- rbind(ptreps, NSdata_centroid)
st_write(spData, "comdata.shp")


# do some data cleaning
tapply(spData$LASTOBS_YR, spData$SUBNATION, summary)  




########################################
# get the study area data ################################################################################################
studyArea <- arc.open(studyArea)
studyArea <- arc.select(studyArea)
studyArea <- arc.data2sf(studyArea)


