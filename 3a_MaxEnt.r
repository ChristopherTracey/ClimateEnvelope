
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

