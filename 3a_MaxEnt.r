####################

# split into groups for k-fold cross validation
group <- kfold(species_pts,5)
pres_train <- species_pts[group!=1,]
pres_test <- species_pts[group==1,]

# MaxEnt Model ########################################################################################################

ifelse(!dir.exists(here::here("_data","species",sp_code,"output","maxent")), dir.create(here::here("_data","species",sp_code,"output","maxent")), FALSE)

outPath <- here::here("_data","species",sp_code,"output","maxent")

me.out <- maxent(predictors_Current, pres_train, path=outPath)

# remove the least important variables
me.dat <- as.data.frame(slot(me.out, "results"))
me.imp.dat <- me.dat[grepl("permutation.importance",rownames(me.dat)), ,drop = FALSE]
me.imp.dat <- cbind(me.imp.dat, "var" = unlist(lapply(rownames(me.imp.dat), FUN = function(x) strsplit(x, "\\.")[[1]][[1]])))
impvals <- me.imp.dat
names(impvals) <- c("imp","var")

OriginalNumberOfEnvars <- nrow(impvals)

db_cem <- dbConnect(SQLite(), dbname=nm_db_file) # connect to the database
corrdEVs <- dbGetQuery(db_cem, "SELECT rasCode, CorrGroup FROM lkpEnvVar WHERE CorrGroup  IS NOT NULL order by CorrGroup ;" )
dbDisconnect(db_cem)

corrdEVs <- corrdEVs[corrdEVs$rasCode %in% impvals$var,]
if(nrow(corrdEVs) > 0 ){
  for(grp in unique(corrdEVs$CorrGroup)){
    vars <- corrdEVs[corrdEVs$CorrGroup == grp,"rasCode"]
    imp.sub <- impvals[impvals$var %in% vars,, drop = FALSE]
    suppressWarnings(varsToDrop <- imp.sub[!imp.sub$imp == max(imp.sub$imp),, drop = FALSE])
    impvals <- impvals[!impvals$var %in% varsToDrop$var,,drop = FALSE]
  }
  rm(vars, imp.sub, varsToDrop)
}











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
)

p_thresh <- reclassify(p, m)
pfut_thresh <- reclassify(pfut, m1)

# creates a map showing contraction/stability/expansion ###################
m3 <- calc(stack(p_thresh,pfut_thresh), sum)
plot(m3)

library(rasterVis)
gplot(m3) +
  geom_raster(aes(fill = factor(value))) +
  coord_equal()
###########################################################################


# p_thresh <- p[p > tr]

plot(p_thresh, main='presence/absence')
plot(pfut_thresh)

