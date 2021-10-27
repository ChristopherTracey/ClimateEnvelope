


# create an output directory if it doesn't exist
ifelse(!dir.exists(here::here("_data","species",sp_code)), dir.create(here::here("_data","species",sp_code)), FALSE)

# get species information from the database
db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
sp_data <- dbGetQuery(db_cem, paste("SELECT * FROM lkpSpecies WHERE CUTECODE = " , sQuote(sp_code), sep="") )
dbDisconnect(db_cem)

model_run_name <- paste0(sp_code, "_" , gsub(" ","_",gsub(c("-|:"),"",as.character(Sys.time()))))


# get the predictor data ##############################################################################################
cat("Loading the predictor data...")

# get the study area
studyArea <- arc.open(studyArea)
studyArea <- arc.select(studyArea)
studyArea <- arc.data2sf(studyArea)

# stack the predictors
predictors_Current <- stack(list.files(pathPredictorsCurrent, pattern = 'tif$', full.names=TRUE ))
predictors_Current <- projectRaster(predictors_Current, crs=crs(studyArea))
predictors_Future <- stack(list.files(pathPredictorsFuture, pattern = 'tif$', full.names=TRUE ))
predictors_Future <- projectRaster(predictors_Future, crs=crs(studyArea))

# check to see if the names are the same
setdiff(names(predictors_Current), names(predictors_Future))
setdiff(names(predictors_Future), names(predictors_Current))

# get set of random points for correlation analysis
set.seed(25)
bg <- dismo::randomPoints(predictors_Current, 10000)
bg <- prepareSWD(species="Bgs", a=bg, env=predictors_Current)
plotCor(bg, method="spearman", cor_th=0.7)
corVar(bg, method="spearman", cor_th=0.7)

# get the species data ################################################################################################
spData <- arc.open(spData_path)
spData <- arc.select(spData) #, where_clause=paste("SNAME = ", sp_data$SNAME, sep="")
spData <- arc.data2sf(spData)

spData_pro <- st_transform(spData, crs=crs(predictors_Current))

# write a shapeefile of the training data to the input folder for backup or other use
ifelse(!dir.exists(here::here("_data","species",sp_code,"input")), dir.create(here::here("_data","species",sp_code,"input")), FALSE)
ifelse(!dir.exists(here::here("_data","species",sp_code,"output")), dir.create(here::here("_data","species",sp_code,"output")), FALSE)
# write.csv(species_pts, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.csv")))
st_write(spData, here::here("_data", "species", sp_code, "input", paste0(sp_code,"_input.shp")), append=FALSE)

# get some metadata
md_ptTraining <- nrow(spData)

# convert species points to a lat/long df  MOVE THIS BELOW???
coords_pres <- data.frame(st_coordinates(spData_pro[,1]))
names(coords_pres) <- c("x","y")

# background coordinates
# These are 500 random locations, used as in place of absence values as 
# 'pseudoabsences' (the species probably doesn't occur at any random point)
coords_bg <- as.data.frame(dismo::randomPoints(predictors_Current, 500))

# # check to see if the coordinate systems are the same
# a <- stringr::str_split(as.character(crs(predictors_Current)), ' ')
# b <- stringr::str_split(as.character(crs(spData)), ' ')
# identical(sort(unlist(a)), sort(unlist(b)))
# rm(a,b)

# create SDW object ############################
data.SWD <- prepareSWD(species=spData$SNAME, p=coords_pres, a=coords_bg, env=predictors_Current)
data.SWD

i = 2
# run the models ##################################

#for(i in 1:length(ModelMethods)){
cat(paste("Running a ", ModelMethods[i], " model for ", unique(spData$SNAME),".", sep=""))
model <- train(method=ModelMethods[i], data=data.SWD)  # model <- train(method = "Maxent", data = data, fc = "lh", reg = 0.5, iter = 700)
model
map <- predict(model, data=predictors_Current, type="cloglog")
plotPred(map)

cat("- variable importance")
if(ModelMethods[i]=="Maxent"){
  vi <- maxentVarImp(model)
  vi
  plotVarImp(vi[, 1:2])
} else if(ModelMethods[i]=="RF"|ModelMethods[i]=="BRT") {
  vi <- SDMtune::varImp(model) 
  vi
  plotVarImp(vi[, 1:2])   
} else {
  cat("No valid variable importance method exists...")
}

# set up folds for training and test
folds <- randomFolds(data.SWD, k=4, only_presence=TRUE, seed=5)
cv_model <- train(ModelMethods[i], data=data.SWD, folds=folds)
cv_model
auc(cv_model)
auc(cv_model, test=TRUE)

# variable selection
selected_variables_model <- varSel(cv_model, metric="auc", bg4cor=bg, method="spearman", cor_th=0.7, permut=1)
selected_variables_model

varImp(model, permut=1)
cat("Testing AUC before: ", auc(cv_model, test=TRUE))
reduced_variables_model <- reduceVar(cv_model, th=6, metric="auc", permut=1)
cat("Testing AUC after: ", auc(reduced_variables_model, test = TRUE))

# jacknifing approach to variable selection
cat("Testing AUC before: ", tss(cv_model, test=TRUE))
reduced_variables_model <- reduceVar(cv_model, th=15, metric="tss", permut=1, use_jk=TRUE)
cat("Testing AUC after: ", tss(reduced_variables_model, test=TRUE))

#}



# run a model #############################################################
model <- train(method="Maxent", data=data)
model

slotNames(model@model)

model <- train(method = "Maxent", data = data, fc = "lh", reg = 0.5, iter = 700)
model

pred <- predict(model, data=data, type="cloglog")
head(pred)

map <- predict(model, data=predictors_Current, type="cloglog")
plotPred(map)


auc(model)
plotROC(model)
tss(model)
aicc(model, env=predictors_Current)

# library(zeallot)  # For unpacking assignment
# c(train, test) %<-% trainValTest(data, test = 0.2, only_presence = TRUE,
#                                  seed = 25)


# model <- train("Maxent", data = train)
# auc(model)
# auc(model, test = test)
# plotROC(model, test = test)

folds <- randomFolds(data, k = 4, only_presence = TRUE, seed = 25)
cv_model <- train("Maxent", data = data, folds = folds)
cv_model
auc(cv_model)
auc(cv_model, test = TRUE)



model@model@results


vi <- maxentVarImp(model)
vi
plotVarImp(vi[, 1:2])

# variable selection

set.seed(25)
bg <- dismo::randomPoints(predictors_Current, 10000)
bg <- prepareSWD(species = "Bgs", a=bg, env=predictors_Current)

plotCor(bg, method="spearman", cor_th=0.7)
corVar(bg, method="spearman", cor_th=0.7)

selected_variables_model <- varSel(model, metric="auc", test=test, bg4cor=bg, method="spearman", cor_th=0.7, permut=1)

selected_variables_model

varImp(model, permut = 1)
cat("Testing AUC before: ", auc(model, test = test))
reduced_variables_model <- reduceVar(model, th = 6, metric = "auc", test=test, permut = 1)
cat("Testing AUC after: ", auc(reduced_variables_model, test = test))


cat("Testing AUC before: ", auc(model, test = test))
reduced_variables_model <- reduceVar(model, th = 15, metric = "auc", test = test, permut = 1, use_jk = TRUE)
cat("Testing AUC after: ", auc(reduced_variables_model, test = test))










###########################################


# run a model #############################################################
model <- train(method="BRT", data=data)
model

slotNames(model@model)

model <- train(method = "BRT", data = data, fc = "lh", reg = 0.5, iter = 700)
model

pred <- predict(model, data=data, type="cloglog")
head(pred)

map <- predict(model, data=predictors_Current, type="cloglog")
plotPred(map)

map1 <- predict(model, data=predictors_Future, type="cloglog")
plotPred(map1)

auc(model)
plotROC(model)
tss(model)

aicc(model, env=predictors_Current)

# library(zeallot)  # For unpacking assignment
# c(train, test) %<-% trainValTest(data, test = 0.2, only_presence = TRUE,
#                                  seed = 25)


# model <- train("Maxent", data = train)
# auc(model)
# auc(model, test = test)
# plotROC(model, test = test)

folds <- randomFolds(data, k = 4, only_presence = TRUE, seed = 25)
cv_model <- train("BRT", data = data, folds = folds)
cv_model
auc(cv_model)
auc(cv_model, test = TRUE)



model@model@results


vi <- SDMtune::varImp(model) # maxentVarImp(model)
vi
plotVarImp(vi[, 1:2])

# variable selection

set.seed(25)
bg <- dismo::randomPoints(predictors_Current, 10000)
bg <- prepareSWD(species = "Bgs", a=bg, env=predictors_Current)

plotCor(bg, method="spearman", cor_th=0.7)
corVar(bg, method="spearman", cor_th=0.7)

selected_variables_model <- varSel(model, metric="auc", test=test, bg4cor=bg, method="spearman", cor_th=0.7, permut=1)

selected_variables_model

varImp(model, permut = 1)
cat("Testing AUC before: ", auc(model, test = test))
reduced_variables_model <- reduceVar(model, th = 6, metric = "auc", test=test, permut = 1)
cat("Testing AUC after: ", auc(reduced_variables_model, test = test))


cat("Testing AUC before: ", auc(model, test = test))
reduced_variables_model <- reduceVar(model, th = 15, metric = "auc", test = test, permut = 1, use_jk = TRUE)
cat("Testing AUC after: ", auc(reduced_variables_model, test = test))



############################## RF



