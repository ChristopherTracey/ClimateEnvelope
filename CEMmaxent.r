library(here)
library(dismo)
library(arcgisbinding)
arc.check_product()
library(ggplot2)

# get the basemap data ################################################################################################
bndPA <- arc.open(here::here("modeling_data.gdb", "bnd_PAstate"))
bndPA <- arc.select(bndPA)
bndPA <- arc.data2sf(bndPA)

# get the predictor data ##############################################################################################
predictors_Current <- stack(list.files(here::here("Climate_Current"), pattern = 'asc$', full.names=TRUE ))
predictors_Future <- stack(list.files(here::here("Climate_Future2050rfp45"), pattern = 'asc$', full.names=TRUE ))
# check to see if the names are the sames
setdiff(names(predictors_Current), names(predictors_Future))
setdiff(names(predictors_Future), names(predictors_Current))

# set the modeling extent to be the extent of the rasters
ext = extent(predictors_Current)

# get the training data ###############################################################################################
species <- read.table(here::here("Maxent_practiceinputspecies.csv"), header=TRUE,  sep=',')
unique(species$species)
species <- species[which(species$species=="Lupinus perennis"),]
species_pts <- species[,-1] # drop the 


# MaxEnt Model ########################################################################################################
# split into groups for k-fold cross validation
group <- kfold(species_pts,5)
pres_train <- species_pts[group!=1,]
pres_test <- species_pts[group==1,]

# make the model
xm <- maxent(predictors_Current, pres_train)
plot(xm) # plot variable importance

response(xm) # plot the response curves

# make random background points
backg = randomPoints(predictors_Current, n=1000, p=species_pts, ext=ext, extf=1.25, excludep=TRUE) # it appears we can add a bias grid into this...
colnames(backg) <- c("lon","lat")
group=kfold(backg, 5)
backg_train <- backg[group!=1,]
backg_test <- backg[group==1,]

e = evaluate(pres_test, backg_test, xm, predictors_Current)
plot(e,'ROC')

# predict to current
p <- predict(predictors_Current, xm, ext=ext, progress='')
# plot(p, main="Maxent, raw values - Present")
# plot(sf::st_geometry(bndPA), add=TRUE)





# predict to the future!
pfut <- predict(predictors_Future, xm, ext=ext, progress='')
# plot(pfut, main="Maxent, raw values - Future")
# plot(sf::st_geometry(bndPA), add=TRUE)

# make a graph
par(mfrow=c(1,2))
plot(p, main='MaxEnt, Current')
plot(sf::st_geometry(bndPA), add=TRUE)
plot(pfut, main='MaxEnt, Future')
plot(sf::st_geometry(bndPA), add=TRUE)
par(mfrow = c(1,1)) # reset to 1x1 plot






tr <- threshold(e, 'spec_sens')

m <- cbind(
  from = c(-Inf, tr),
  to = c(tr, Inf),
  becomes = c(0, 1)
)

m1 <-  cbind(
  from = c(-Inf, tr),
  to = c(tr, Inf),
  becomes = c(0, 2)
)µ

p_thresh <- reclassify(p, m)
pfut_thresh <- reclassify(pfut, m1)


m3 <- calc(stack(p_thresh,pfut_thresh), sum)
plot(m3)

library(rasterVis)
gplot(m3) +
  geom_raster(aes(fill = factor(value))) +
  coord_equal()

# p_thresh <- p[p > tr]

plot(p_thresh, main='presence/absence')
plot(pfut_thresh)


#######################
#RF
train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- extract(predictors_Current, train)
envtrain <- data.frame(cbind(pa=pb_train, envtrain) )
#envtrain[,'biome'] = factor(envtrain[,'biome'], levels=1:14)
head(envtrain)

testpres <- data.frame( extract(predictors_Current, pres_test) )
testbackg <- data.frame( extract(predictors_Current, backg_test) )

# random forest
library(randomForest)
model <- pa ~ CMD + DD5 + MAP + MAR + MSP + PAS + RH + TD
rf1 <- randomForest(model, data=envtrain, na.action=na.roughfix)
model <- factor(pa) ~ CMD + DD5 + MAP + MAR + MSP + PAS + RH + TD
rf2 <- randomForest(model, data=envtrain, na.action=na.roughfix)
rf3 <- randomForest(envtrain[,1:8], factor(pb_train))
erf <- evaluate(testpres, testbackg, rf1)
erf

pr <- predict(predictors_Current, rf1, ext=ext)
prf <- predict(predictors_Future, rf1, ext=ext)
par(mfrow=c(1,2))
plot(pr, main='Random Forest, Current')
plot(sf::st_geometry(bndPA), add=TRUE)
plot(prf, main='Random Forest, Future')
plot(sf::st_geometry(bndPA), add=TRUE)
par(mfrow = c(1,1)) # reset to 1x1 plot

# plot(wrld_simpl, add=TRUE, border='dark grey')
# tr <- threshold(erf, 'spec_sens')
# plot(pr > tr, main='presence/absence')
# plot(wrld_simpl, add=TRUE, border='dark grey')
# points(pres_train, pch='+')
# points(backg_train, pch='-', cex=0.25)



