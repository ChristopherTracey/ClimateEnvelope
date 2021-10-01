
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



