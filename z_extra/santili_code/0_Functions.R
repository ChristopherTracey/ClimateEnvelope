##########################################################################################################################################
#                                                                                                                                        #
# Author: XXXX XXXX                                                                                                                	     #
# e-mail: XXXXXX@XXXXX                                                                                                    				 #
# Reference:                                                                                                                             #
# XXXX et al. 2019. Assessing the reliability of species distribution projections in climate change research. XXXXXX                     #
#                                                                                                                                        #
##########################################################################################################################################

#This code includes a collection of functions used by "Run_SDMs.R"

require(pROC)
require(PresenceAbsence)
require(maxnet)
require(randomForest)
require(dismo)
require(rgeos)
require(mgcv)
require(gbm)
require(caret)
require(raster)
require(hypervolume)
require(virtualspecies)
require(usdm)
require(GISTools)
require(MASS)

vifstep2<-function(x, tr=3){
excluded<-c()
v<-vif(x)
while(max(v$VIF)>tr){
varExcl<-as.character(v[v$VIF==max(v$VIF), 'Variables'])
excluded<-c(excluded,varExcl)
x<-x[,which(!names(x) %in% varExcl)]
v<-vif(x)
}
return(excluded)
}

#MODIFIED VERSION OF FORMATFUNCTION THAT TAKES A LIST INSTEAD OF VARIABLES DIRECTLY
formatFunctions_list <- function(x = NULL, rescale = TRUE, listArgs)
{
  details <- list()
  args <- listArgs
  for (i in names(args))
  {
    if(!("fun" %in% names(args[[i]])))
    {
      stop(paste("No function was correctly provided for variable", i))
    }
    details[[i]]$fun <- args[[i]]["fun"]
    args[[i]] <- args[[i]][-(names(args[[i]]) %in% "fun")]
    details[[i]]$args <- as.list(args[[i]])
    a <- sapply(args[[i]], as.numeric)
    if(any(is.na(a)))
    {
      details[[i]]$args[[which(!is.na(a))]] <- sapply(details[[i]]$args[[which(!is.na(a))]], as.numeric)
    } else
    {
      details[[i]]$args <- a
    }
  }
  if (is(x, "Raster"))
  {
    plotResponse(x = x, parameters = details, rescale = rescale, approach = "response")
  } else if (!is.null(x))
  {
    stop("x must either be NULL or a raster stack of environmental variables")
  }
  return(details)
}


#FUNCTION TO GENERATE A RANDOM SPECIES IN A RANDOM CONTINENT

#n = number of points to sample to estimate the niche
#minCell = min number of land cells in the extent sampled
#LandBlocks = raster of land land blocks
#s = raster stack with env variables
#RandomExtent = if true the extent is chosen randomly
#SamplingExtent = if random extent is false, this indicates the size of each side of the extent from which to sample


GenerateRandomSP_19bioclim=function(LandBlocks=LandBlocks, nImpVar=6, RandomExtent=TRUE, SamplingExtent=5, stack_pres=s, stack_future=s2, stack_PCA=sPCA, n=100, minCell=5, randomContinent=TRUE, CONT=87, PAlogistic=FALSE) {

if(randomContinent==TRUE) {
CONT<-sample(c(EU, AF, AU, NAm, SA), 1)
}
#only retain important variables
ImpVar<-sample(names(stack_pres), nImpVar)
stack_pres<-stack_pres[[ImpVar]]
stack_future<-stack_future[[ImpVar]]

#select land block
LB=LandBlocks
LB[LB!=CONT]<-NA
LB<-trim(LB) #crop to non-NA values

stack_pres<-crop(stack_pres, LB)
stack_pres<-mask(stack_pres, LB)
stack_future<-crop(stack_future, LB)
stack_future<-mask(stack_future, LB)

repeat { #repeat sampling of extent until meets conditions

pp<-rasterToPoints(LB)
pp<-pp[sample(nrow(pp), 1), c(1,2)] #random point

if(RandomExtent==TRUE) {
	mn=3#deg
	mx=10
	SamplingExtent<-runif(1, mn, mx)
	xmin=pp[1]-(SamplingExtent/2)
	xmax=pp[1]+(SamplingExtent/2)
	ymin=pp[2]-(SamplingExtent/2)
	ymax=pp[2]+(SamplingExtent/2)

#limit to continent
	xmin<-ifelse(xmin<extent(LB)[1], extent(LB)[1], xmin)
	xmax<-ifelse(xmax>extent(LB)[2], extent(LB)[2], xmax)
	ymin<-ifelse(ymin<extent(LB)[3], extent(LB)[3], ymin)
	ymax<-ifelse(ymax>extent(LB)[4], extent(LB)[4], ymax)
} else { #uses SamplingExtent to extend the extent around the random point
xmin=pp[1]-(SamplingExtent/2)
xmax=pp[1]+(SamplingExtent/2)
ymin=pp[2]-(SamplingExtent/2)
ymax=pp[2]+(SamplingExtent/2)

#limit to continent
	xmin<-ifelse(xmin<extent(LB)[1], extent(LB)[1], xmin)
	xmax<-ifelse(xmax>extent(LB)[2], extent(LB)[2], xmax)
	ymin<-ifelse(ymin<extent(LB)[3], extent(LB)[3], ymin)
	ymax<-ifelse(ymax>extent(LB)[4], extent(LB)[4], ymax)
}
StudyArea<-c(xmin, xmax, ymin, ymax)
#plot(LB, box=FALSE, axes=FALSE); abline(v=xmin); abline(v=xmax); abline(h=ymin); abline(h=ymax)
LBcrop<-crop(LB, StudyArea)
fr<-as.data.frame(freq(LBcrop))
fr<-fr[!is.na(fr$value),]
if(nrow(fr)>0) {if(fr$count>minCell) {break}} #does it include at least n cells of land?
}
#n cells in the study area divided by X (so density of points is constant)
n<-ceiling(fr$count/10)
pp<-sampleRandom(stack_pres, n, na.rm=TRUE, ext=StudyArea, xy=TRUE) #sample random points within the extent
#points(pp[,1:2], cex=0.1, pch=19)

mn_values<-apply(pp[,3:ncol(pp)], 2, mean)
sd_values<-apply(pp[,3:ncol(pp)], 2, sd)

#If not important the SD is replaced with a huge sd so the relationship becomes almost flat

bio_01<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_01']), sd = as.numeric(sd_values['bio_01']))
bio_02<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_02']), sd = as.numeric(sd_values['bio_02']))
bio_03<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_03']), sd = as.numeric(sd_values['bio_03']))
bio_04<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_04']), sd = as.numeric(sd_values['bio_04']))
bio_05<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_05']), sd = as.numeric(sd_values['bio_05']))
bio_06<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_06']), sd = as.numeric(sd_values['bio_06']))
bio_07<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_07']), sd = as.numeric(sd_values['bio_07']))
bio_08<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_08']), sd = as.numeric(sd_values['bio_08']))
bio_09<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_09']), sd = as.numeric(sd_values['bio_09']))
bio_10<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_10']), sd = as.numeric(sd_values['bio_10']))
bio_11<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_11']), sd = as.numeric(sd_values['bio_11']))
bio_12<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_12']), sd = as.numeric(sd_values['bio_12']))
bio_13<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_13']), sd = as.numeric(sd_values['bio_13']))
bio_14<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_14']), sd = as.numeric(sd_values['bio_14']))
bio_15<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_15']), sd = as.numeric(sd_values['bio_15']))
bio_16<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_16']), sd = as.numeric(sd_values['bio_16']))
bio_17<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_17']), sd = as.numeric(sd_values['bio_17']))
bio_18<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_18']), sd = as.numeric(sd_values['bio_18']))
bio_19<-c(fun = 'dnorm', mean = as.numeric(mn_values['bio_19']), sd = as.numeric(sd_values['bio_19']))

LIST_VAR<-list()
for (i in 1:length(ImpVar)){LIST_VAR[[i]]<-get(ImpVar[i])}; names(LIST_VAR)<-ImpVar
my.parameters <- formatFunctions_list(NULL, listArgs=LIST_VAR)

#predict present
VS <- suppressWarnings(generateSpFromFun(raster.stack = stack_pres, parameters = my.parameters, plot=FALSE, species.type='additive'))
#plotResponse(VS)
suit<-VS$suitab.raster #extract raster
#plot(suit, box=FALSE, axes=FALSE)
suit_cr<-crop(suit, StudyArea)

if(PAlogistic){
qq<-quantile(suit, c(0.2, 0.8))
beta<-runif(1, qq[1], qq[2])
PA <- suppressWarnings(convertToPA(VS, beta = beta, alpha = -0.05, plot=FALSE)) #beta and alpha determine the shape of relationship between suitability and probability of presence
PA<-PA$pa.raster #extract raster
#plot(PA, box=FALSE, axes=FALSE)
#CROP PRESENCES (BUT NOT SUITABILITY) TO THE ORIGINAL STUDY AREA
PA<-crop(PA, StudyArea)
#plot(PA, box=FALSE, axes=FALSE)
} else {
qq<-quantile(suit_cr, c(0.3, 0.9)) #from 10 a 70% of prevalence
tr<-runif(1, qq[1], qq[2])
PA<-suit
PA[PA>tr]<-1
PA[PA<tr]<-0
PAcr<-crop(PA, StudyArea) #PA has potential presences across the continent, PAcr has the realized ones
}


#predict future
VS2 <- generateSpFromFun(raster.stack = stack_future, parameters = my.parameters, plot=FALSE, species.type='additive')
suit2<-VS2$suitab.raster #extract raster
suit2_cr<-crop(suit2, StudyArea)

#plot(suit2, box=FALSE, axes=FALSE)
if(PAlogistic){
PA2 <- convertToPA(VS2, beta = beta, alpha = -0.05, plot=FALSE) #beta and alpha determine the shape of relationship between suitability and probability of presence
PA2<-PA2$pa.raster #extract raster
#plot(PA2, box=FALSE, axes=FALSE)
#CROP FUTURE PRESENCES (BUT NOT SUITABILITY) TO THE ORIGINAL STUDY AREA
PA2<-crop(PA2, StudyArea)
#plot(PA2, box=FALSE, axes=FALSE)
} else {
PA2<-suit2
PA2[PA2>tr]<-1
PA2[PA2<tr]<-0
PA2cr<-crop(PA2, StudyArea)
}



#Estimate changes
delta<-PA2-PA
FREQd<-as.data.frame(freq(delta))
FREQd<-FREQd[!is.na(FREQd$value),]
loss<-FREQd[FREQd$value==-1,'count']
gain<-FREQd[FREQd$value==1,'count']

FREQpa<-as.data.frame(freq(PA))
FREQpa<-FREQpa[!is.na(FREQpa$value),]
AreaPresent<-FREQpa[FREQpa$value==1,'count']

FREQpa2<-as.data.frame(freq(PA2))
FREQpa2<-FREQpa2[!is.na(FREQpa2$value),]
AreaFuture<-FREQpa2[FREQpa2$value==1,'count']

gain_per<-100*gain/AreaPresent
loss_per<-100*loss/AreaPresent

stats=data.frame(Var=c('AreaPresent', 'AreaFuture', 'gain', 'gain_per', 'loss', 'loss_per', 'HyperVolume', 'SpeciesPrevalence'), Value=c(AreaPresent, AreaFuture, gain, gain_per, loss, loss_per, HyperVolume, tr))

OUT<-list(suit_pt=suit, PA_pt=PA, suit_rl=suit_cr, PA_rl=PAcr, suit2_pt=suit2, PA2_pt=PA2, suit2_rl=suit2_cr, PA2_rl=PA2cr, LB=LB, stats=stats, ImpVar=ImpVar)
return(OUT)
}

GenerateRandomSP=function(LandBlocks=LandBlocks, RandomExtent=TRUE, SamplingExtent=5, stack_pres=s, stack_future=s2, n=100, minCell=5, randomContinent=TRUE, CONT=87, HypervolumeGaussian=FALSE, Plot=TRUE) {

if(randomContinent==TRUE) {
CONT<-sample(c(EU, AF, AU, NAm, SA), 1)
}

LB=LandBlocks
LB[LB!=CONT]<-NA
LB<-trim(LB) #crop to non-NA values

stack_pres<-crop(stack_pres, LB)
stack_pres<-mask(stack_pres, LB)
stack_future<-crop(stack_future, LB)
stack_future<-mask(stack_future, LB)

#Create species
#set extent
repeat { #repeat sampling of extent until meets conditions

pp<-rasterToPoints(LB)
pp<-pp[sample(nrow(pp), 1), c(1,2)] #random point

if(RandomExtent==TRUE) {
	mn=3#deg
	mx=50
	SamplingExtent<-runif(1, mn, mx)
	xmin=pp[1]-(SamplingExtent/2)
	xmax=pp[1]+(SamplingExtent/2)
	ymin=pp[2]-(SamplingExtent/2)
	ymax=pp[2]+(SamplingExtent/2)

#limit to continent
	xmin<-ifelse(xmin<extent(LB)[1], extent(LB)[1], xmin)
	xmax<-ifelse(xmax>extent(LB)[2], extent(LB)[2], xmax)
	ymin<-ifelse(ymin<extent(LB)[3], extent(LB)[3], ymin)
	ymax<-ifelse(ymax>extent(LB)[4], extent(LB)[4], ymax)
#xmin=runif(1, extent(LB)[1], pp[1])
#xmax=runif(1, pp[1], extent(LB)[2])
#ymin=runif(1, extent(LB)[3], pp[2])
#ymax=runif(1, pp[2], extent(LB)[4])
} else { #uses SamplingExtent to extend the extent around the random point
xmin=pp[1]-(SamplingExtent/2)
xmax=pp[1]+(SamplingExtent/2)
ymin=pp[2]-(SamplingExtent/2)
ymax=pp[2]+(SamplingExtent/2)

#limit to continent
	xmin<-ifelse(xmin<extent(LB)[1], extent(LB)[1], xmin)
	xmax<-ifelse(xmax>extent(LB)[2], extent(LB)[2], xmax)
	ymin<-ifelse(ymin<extent(LB)[3], extent(LB)[3], ymin)
	ymax<-ifelse(ymax>extent(LB)[4], extent(LB)[4], ymax)
}
StudyArea<-c(xmin, xmax, ymin, ymax)
#plot(LB); abline(v=xmin); abline(v=xmax); abline(h=ymin); abline(h=ymax)
LBcrop<-crop(LB, StudyArea)
fr<-as.data.frame(freq(LBcrop))
fr<-fr[!is.na(fr$value),]
if(nrow(fr)>0) {if(fr$count>minCell) {break}} #does it include at least n cells of land?
}
#n cells in the study area divided by X (so density of points is constant)
n<-ceiling(fr$count/10)
pp<-sampleRandom(stack_pres, n, na.rm=TRUE, ext=StudyArea, xy=TRUE) #sample random points within the extent

#estimate mean and sd of the area for each variable (otherwise is difficult to come up with something realistic)
mn_values<-apply(pp[,3:4], 2, mean)
sd_values<-apply(pp[,3:4], 2, sd)

#set niche parameters 4 bioclim variables
#my.parameters <- formatFunctions(Tm = c(fun = 'dnorm', mean = as.numeric(mn_values['Tm']), sd = as.numeric(sd_values['Tm'])),
#                                 Tsd = c(fun = 'dnorm', mean = as.numeric(mn_values['Tsd']), sd = as.numeric(sd_values['Tsd'])),
#                                 Pm = c(fun = 'dnorm', mean = as.numeric(mn_values['Pm']), sd = as.numeric(sd_values['Pm'])),
#                                 Pcv = c(fun = 'dnorm', mean = as.numeric(mn_values['Pcv']), sd = as.numeric(sd_values['Pcv']))
#                                 )

#set niche parameters 2 PCA axes
my.parameters <- formatFunctions(PC1 = c(fun = 'dnorm', mean = as.numeric(mn_values['PC1']), sd = as.numeric(sd_values['PC1'])),
                                 PC2 = c(fun = 'dnorm', mean = as.numeric(mn_values['PC2']), sd = as.numeric(sd_values['PC2']))
                                 )

#stack_pres<-crop(stack_pres, StudyArea)
#stack_future<-crop(stack_future, StudyArea)

#predict present
VS <- suppressWarnings(generateSpFromFun(raster.stack = stack_pres, parameters = my.parameters, plot=FALSE))
#plotResponse(VS)
suit<-VS$suitab.raster #extract raster

PA <- suppressWarnings(convertToPA(VS, beta = 0.5, alpha = -0.05, plot=FALSE)) #beta and alpha determine the shape of relationship between suitability and probability of presence
PA<-PA$pa.raster #extract raster

#CROP PRESENCES (BUT NOT SUITABILITY) TO THE ORIGINAL STUDY AREA
tmp<-crop(PA, StudyArea)
tmp<-extend(tmp, PA)
tmp[is.na(tmp)]<-0
PA[!is.na(PA)]<-0
PA<-PA+tmp
#crop it to the continent
#suit<-crop(suit, extent(LB)) 
#suit<-mask(suit, LB)

#PA<-crop(PA, extent(LB)) 
#PA<-mask(PA, LB)

#predict future
VS2 <- generateSpFromFun(raster.stack = stack_future, parameters = my.parameters, plot=FALSE)
suit2<-VS2$suitab.raster #extract raster

PA2 <- convertToPA(VS2, beta = 0.5, alpha = -0.05, plot=FALSE) #beta and alpha determine the shape of relationship between suitability and probability of presence
PA2<-PA2$pa.raster #extract raster

#CROP FUTURE PRESENCES (BUT NOT SUITABILITY) TO THE ORIGINAL STUDY AREA
tmp<-crop(PA2, StudyArea)
tmp<-extend(tmp, PA2)
tmp[is.na(tmp)]<-0
PA2[!is.na(PA2)]<-0
PA2<-PA2+tmp

#crop it to the continent
#suit2<-crop(suit2, extent(LB)) 
#suit2<-mask(suit2, LB)

#PA2<-crop(PA2, extent(LB)) 
#PA2<-mask(PA2, LB)

#Hypervolume of the niche
PAmask<-PA
PAmask[PAmask!=1]<-NA
stack_pres_scale<-scale(stack_pres)
stack_cropped<-crop(stack_pres_scale, extent(PAmask))
EnvInPA<-mask(stack_cropped, PAmask)
EnvInPA<-rasterToPoints(EnvInPA)[,-c(1:2)]

if(HypervolumeGaussian==TRUE) {hv<-hypervolume_gaussian(EnvInPA, verbose=FALSE)} else {hv<-hypervolume_box(EnvInPA, verbose=FALSE)}
HyperVolume<-hv@Volume

if(Plot==TRUE) {
par(mfrow=c(2,2))
plot(suit, main='Present suitability')
plot(PA, main='Present suitability')
plot(suit2, main='Future suitability')
plot(PA2, main='Future Presence')
}

#Estimate changes
delta<-PA2-PA
loss<-freq(delta)[1,'count']
gain<-freq(delta)[3,'count']
AreaPresent<-freq(PA)[2,'count']
AreaFuture<-freq(PA2)[2,'count']
gain_per<-100*gain/AreaPresent
loss_per<-100*loss/AreaPresent

stats=data.frame(Var=c('AreaPresent', 'AreaFuture', 'gain', 'gain_per', 'loss', 'loss_per', 'HyperVolume'), Value=c(AreaPresent, AreaFuture, gain, gain_per, loss, loss_per, HyperVolume))

OUT<-list(suit=suit, PA=PA, suit2=suit2, PA2=PA2, LB=LB, stats=stats)
return(OUT)
}

#split dataset in training and testing
SplitDataset<-function(data, partition) {

Pres<-data[data$PA==1,]
Ab<-data[data$PA==0,]
nrPres<-nrow(Pres)
nrAb<-nrow(Ab)
trainRowsP<-sample(nrPres, (partition/100)*nrPres) #sample n random rows of the dataset
trainRowsA<-sample(nrAb, (partition/100)*nrAb) #sample n random rows of the dataset
trainingP=Pres[trainRowsP,] #extract the rows to generate the training dataset
testingP=Pres[-trainRowsP,] #use the remnant rows to generate the testing dataset
trainingA=Ab[trainRowsA,] #extract the rows to generate the training dataset
testingA=Ab[-trainRowsA,] #use the remnant rows to generate the testing dataset
training=rbind(trainingP, trainingA)
testing=rbind(testingP, testingA)
return(list(training, testing))
}

#assign percentage area to the polygons of a multipolygon
AreaPercent<-function(x) {
  tot_area <- sum(sapply(slot(x, "polygons"), slot, "area"))
  x2<-sapply(slot(x, "polygons"), slot, "area") / tot_area * 100
  return(x2)
}

PercentageBuffer <- function(f, percent, precision, maxsteps = 20){
    require(rgeos)
    A0 = gArea(f)

    targetA = A0 + A0 * (percent/100)
    e = bbox(f)
    w0 = 0
    w2 = diff(e[1,])/10 # initial w

    ## initial expansion until we pass the solution:

    repeat{
        b1 = gBuffer(f, width=w2)
        if(gArea(b1) < targetA){
            w2 = w2 * 2
        }else{
            break
        }
    }
    w1 = (w0 + w2)/2

    ## solution is between w0 and w2

    steps = 1

    repeat{
        b = gBuffer(f, width=w1)
        A = gArea(b)
        if(abs(100*(A-targetA)/targetA) < precision){
            return(list(steps = steps, poly = b, target=targetA))
        }
        if(A < targetA){
            w0 = w1
            w1 = (w1+w2)/2
        }
        if(A > targetA){
            w2 = w1
            w1 = (w2+w0)/2
        }
        steps = steps + 1
        if(steps > maxsteps){
            return(NA)
        }
    }
}

#SAMPLE POINTS AND FIT SDM ON VIRTUAL SPECIES
#PA = PresenceAbsence map from which to sample
#PA2 = PresenceAbsence map of the future
#LB = land block raster
#n_pres = number of presence points
#n_back = number of background data
#SpatVal = spatial block validation 
#replicates = n replicates for model fitting and evaluation/treshold estimation
#average the predictions of different replicates and apply the average threshold
#AbsNotInPres = if true, background points are not sampled in presence cells
#backgroundWithinBuffer=study area is defined using a buffer, if true background points are only sampled within the buffer
#buffP = buffer in %
#backgroundInRange = Valid if rangeSDM=TRUE, if true background are sampled within the range of the species, otherwise are sampled within the buffer only
#BiasMap = Use map for biased sampling in species range? 
#nonEquilibrium = assuming species are in non-equilibrium (potentially suitable areas are unoccupied)
#nonEquilibriumHIIvalue = quantile of HII within species potential range above which the species is absent, only works if nonEquilibrium=TRUE


FitSDM<-function(LB=LB, SUIT=ST, SUIT2=ST2, PA=PA, PA2=PA2, Ext_Rls=Ext_Rls, stack_pres=s, stack_future=s2, truePred=varRet, falsePred=varAdd, 
	Model='MaxEnt', MXPredType='cloglog', ExcludeCollinear=TRUE, equalPrevalence=TRUE, n_pres=100, n_back=10000, SpatVal=TRUE, replicates=10, 
	backgroundWithinBuffer=TRUE, bufferP=20, AbsNotInPres=TRUE, 
	rangeSDM=FALSE, backgroundInRange=TRUE, 
	HyperVolumePresAbs=FALSE, HypervolumeGaussian=FALSE,
	Bias=TRUE, BiasQ=100, BiasMap=NULL, nonEquilibrium=FALSE, HIImap=NULL, nonEquilibriumQuantileValue=0.75,
	Plot=FALSE) 
	{

#crop potentially occupied area by realized niche
PRES_mask<-crop(PA, Ext_Rls)
PRES_mask[PRES_mask==0]<-NA
PRES_mask2<-crop(PA2, Ext_Rls)
PRES_mask2[PRES_mask2==0]<-NA

#if nonEquilibrium=TRUE part of the range is unoccupied because of human impact
if(nonEquilibrium==TRUE) {
	PRES_mask_HII<-crop(HIImap, extent(PRES_mask))
	#plot(PRES_mask_HII, box=FALSE, axes=FALSE)
  PRES_mask_HII<-mask(PRES_mask_HII, PRES_mask)
	t<-quantile(PRES_mask_HII, nonEquilibriumQuantileValue) #75ft percentile of HII in species range
	ncells<-nrow(rasterToPoints(PRES_mask_HII))
	PRES_mask_HII[PRES_mask_HII>t]<-NA
	ncells2<-nrow(rasterToPoints(PRES_mask_HII))
  #plot(PRES_mask, box=FALSE, axes=FALSE, col='red')
  #plot(PRES_mask_HII, col='darkgreen', add=TRUE)
	PRES_mask<-mask(PRES_mask, PRES_mask_HII)
  propNonEq=100*ncells2/ncells #percentage of total potential range that is occupied
}
propNonEq<-ifelse(exists('propNonEq'), propNonEq, 0)

set.seed(1)
if(Bias==FALSE) {
	tmp<-rasterToPoints(PRES_mask)
	if(n_pres > nrow(tmp)) {n_pres = nrow(tmp)} #if n of pres > than n cell, then sample max n cell
	pres_sd <- sampleRandom(PRES_mask, n_pres, xy = TRUE, sp=TRUE, na.rm = TRUE)
	pres<-coordinates(pres_sd)
} else {
	
	#with random env var (true if exist, spurious otherwise)
	if(length(truePred)>0){biasVar<-truePred[sample(length(truePred), 1)]} else {biasVar<-falsePred[sample(length(falsePred), 1)]}
	biasVar<-crop(s[[biasVar]], extent(PRES_mask))
	PRES_mask_bias<-crop(biasVar, extent(PRES_mask))
	PRES_mask_bias<-mask(PRES_mask_bias, PRES_mask)
  	#plot(PRES_mask_bias, box=FALSE, axes=FALSE)
  	q<-quantile(PRES_mask_bias, 1-(BiasQ/100))
  	PRES_mask_bias[PRES_mask_bias<q]<-NA
	tmp<-rasterToPoints(PRES_mask_bias)
	if(n_pres > nrow(tmp)) {n_pres = nrow(tmp)} #if n of pres > than n cell, then sample max n cell
	pres_sd <- sampleRandom(PRES_mask_bias, n_pres, xy = TRUE, sp=TRUE, na.rm = TRUE)
	pres<-coordinates(pres_sd)

	#plot(log10(PRES_mask_bias), box=FALSE, axes=FALSE); points(pres, cex=0.4, pch=19) #PLOT sampled points

  kernel<-kde.points(SpatialPoints(pres),0.5) #fit kernel on sampled points
  kernel<-raster(kernel) #rasterize
  kernel<-resample(kernel, PRES_mask) #re-extent to the whole present area
  kernel[is.na(kernel)]<-0 #convert empty areas to zero
  kernel<-mask(kernel, PRES_mask) #mask
  sdBias<-sd(rasterToPoints(kernel)[,3]) #estimate bias as the sd of the kernel within the present area 
  cvBias<-sd(rasterToPoints(kernel)[,3])/mean(rasterToPoints(kernel)[,3])
}
sdBias<-ifelse(exists('sdBias'), sdBias, 0)
cvBias<-ifelse(exists('cvBias'), cvBias, 0)

if(backgroundWithinBuffer==TRUE) { #define study area for sampling pseudo-absences

tmp<-rasterToPoints(PRES_mask)
tmp<-tmp[tmp[,3]==1,1:2]
pres_mcp<-gConvexHull(SpatialPoints(tmp)) #convex hull around known presences (not presence points)
#plot(PA, box=FALSE, axes=FALSE);plot(pres_mcp, add=TRUE)

buffV<-suppressWarnings(PercentageBuffer(pres_mcp, percent=buffP, precision=0.5))$poly
#plot(buffV, box=FALSE, axes=FALSE);plot(PA3, add=TRUE, col='chartreuse4', legend=FALSE);plot(pres_mcp, add=TRUE)
buff<-rasterize(buffV, LB)
buff<-mask(buff, LB)
buff<-trim(buff)
buffForSampling<-buff #if AbsNotInPres=TRUE then it will sample from here
#Calc env dissimilarity (MESS)
s_masked<-crop(sb, extent(buff))
s_masked<-mask(s_masked, buff)
xys <- extract(s_masked, pres)
ms <- mess(s_masked, xys, full=FALSE)
ms[is.infinite(ms)]<-NA

if(AbsNotInPres==TRUE) {
	pres_points_rast<-rasterFromXYZ(cbind(pres, 1))
	pres_points_rast<-extend(pres_points_rast, buffForSampling)
	buffForSampling<-extend(buffForSampling, pres_points_rast)
	buffForSampling<-mask(buffForSampling, pres_points_rast, inverse=TRUE)
}

if(n_back>nrow(rasterToPoints(buffForSampling))) {n_back=nrow(rasterToPoints(buffForSampling))} #if more background than cells, n background==ncell
set.seed(1)
background <- sampleRandom(buffForSampling, n_back, xy = TRUE, sp=TRUE, na.rm = TRUE)
background<-coordinates(background)
	} else { #sample background across the whole buffer including present cells
if(n_back>nrow(rasterToPoints(LB))) {n_back=nrow(rasterToPoints(LB))} #if more background than cells, n background==ncell
set.seed(1)
background <- sampleRandom(LB, n_back, xy = TRUE, sp=TRUE, na.rm = TRUE)
background<-coordinates(background)
}



#species prevalence calculation
#plot(buff); plot(PRES_mask, col='red', add=TRUE)#Plot species prevalence
ppPA<-table(rasterToPoints(PRES_mask)[,3])['1'] #realized presence cells
ppPA2<-table(rasterToPoints(PRES_mask2)[,3])['1'] #realized future cells
studyArea<-table(rasterToPoints(buff)[,3])['1'] #study area cells
SpPrev<-ppPA/studyArea

#Sp prev within the extent sampled for background points
PA_buff<-mask(PA, buffV)
tb<-table(rasterToPoints(PA_buff)[,3])
pa_ratio<-tb['1']/sum(tb)

#dataset
SpData<-as.data.frame(rbind(cbind(pres, 1), cbind(background, 0)))
names(SpData)<-c('x', 'y', 'PA')

Env<-extract(stack_pres, SpData[,c('x','y')])
SpData<-as.data.frame(cbind(SpData, Env))


#EXCLUDE COLLINEAR VARIABLES
if(ExcludeCollinear==TRUE){
excl<-vifstep2(SpData[,4:ncol(SpData)], tr=3)
SpData<-SpData[,which(!names(SpData) %in% excl)]
}

#env variables for predictions 
#crop to the study area
s_crop<-crop(stack_pres, buff)
s_crop<-mask(s_crop, buff)
EnvPres<-as.data.frame(rasterToPoints(s_crop))
EnvPres<-EnvPres[complete.cases(EnvPres),]
s2_crop<-crop(stack_future, buff)
s2_crop<-mask(s2_crop, buff)
EnvFut<-as.data.frame(rasterToPoints(s2_crop))
EnvFut<-EnvFut[complete.cases(EnvFut),]

#dissimilarity pres and fut conditions
ms<-mess(s_crop, EnvFut[,-c(1,2)])
ms<-rasterToPoints(ms)
ms<-ms[complete.cases(ms) & is.finite(ms[,3]),]
PresFutSimilarity_propMoreThan0=nrow(ms[ms[,3]>0,])/nrow(ms)
PresFutSimilarity_Mean=mean(ms[,3], na.rm=TRUE)

#vectors to fill for cross-validatin
tss_tr_cv<-rep(NA, replicates)
eqss_tr_cv<-rep(NA, replicates)

sens_cv<-rep(NA, replicates)
spec_cv<-rep(NA, replicates)
tss_cv<-rep(NA, replicates)
auc_cv<-rep(NA, replicates)

if(SpatVal==TRUE){
	TrainTest<-SpatBlocking(SpData)
	SpData2<-TrainTest[['data']]
	SpData2<-SpData2[!is.na(SpData2$Block),]
	blocks<-TrainTest[['blocks']]
} else {blocks<-replicates}

for (i in 1:length(blocks)) {

if(SpatVal==TRUE){
	B<-blocks[i]
	training<-SpData2[!SpData2$Block %in% B,]
	testing<-SpData2[SpData2$Block %in% B,]
	training<-training[complete.cases(training),]
	testing<-testing[complete.cases(testing),]
} else {
	TrainTest<-SplitDataset(SpData, partition=80) #with partition=100, no testing dataset is created
	training=TrainTest[[1]] 
	testing=TrainTest[[2]]
	training<-training[complete.cases(training),]
	testing<-testing[complete.cases(testing),]
}

#weights
if(equalPrevalence==TRUE) {
tw<-table(training$PA)#n of pres and abs
wp<-(sum(tw)/2)/tw['1'] #weight pres
wb<-(sum(tw)/2)/tw['0'] #weight abs
weights<-ifelse(training$PA==1, wp, wb)
} else {weights=rep(1, nrow(training))}

if(Model == 'GLM') {
 	F<-as.formula(paste0('PA ~ ', paste(paste0('poly(',names(training[,-c(1:3)]),',2)'), collapse='+'))) #first order polynomial
	mod<-suppressWarnings(stepAIC(glm(F, family='binomial', weights=weights, data=training), trace=FALSE))
	prediction<-predict(mod, testing, type='response')
}
if(Model == 'MaxEnt') {
	mod<-maxnet(p=training$PA, data=training[,-c(1:3)], maxnet.formula(p=training$PA, data=training[,-c(1:3)], classes="lq"))
	prediction<-predict(mod, testing, clamp=TRUE, type=MXPredType)[,1]
}

if(Model == 'RF') {
	F<-as.formula(paste0('PA ~ ', paste(names(training[,-c(1:3)]), collapse='+')))
	trainingRF<-training
	trainingRF$PA<-as.factor(trainingRF$PA)

    mod <- randomForest(formula = F,
                        data = trainingRF,
                        ntree = 500,
                        mtry = round(sqrt(length(names(training[,-c(1:3)])))),
                        importance = FALSE,
                        norm.votes = TRUE,
                        strata = factor(c(0,1))
                        )

 	prediction<-predict(mod, testing, type='prob')[,'1']
}

#CROSS VALIDATION
#tresholds
opt<-findOptimum(predProb=prediction, obs=testing[,'PA'], n_tr=100)

tss_tr_cv[i] <- median(opt[opt$TSS==max(opt$TSS),'Treshold']) #using median because there might be multiple tresholds giving the same TSS
eqss_tr_cv[i] <- median(opt[opt$SSAbsDif==min(opt$SSAbsDif),'Treshold'])

#evaluation
PA_tss<-ifelse(prediction > tss_tr_cv[i], 1, 0)
PA_eqSS<-ifelse(prediction > eqss_tr_cv[i], 1, 0)

EV<-modEval(predProb=prediction, predBin=PA_tss, obs=testing[,'PA'])
sens_cv[i]<-EV['Sensitivity','Value']
spec_cv[i]<-EV['Specificity','Value']
tss_cv[i]<-EV['TSS','Value']
auc_cv[i]<-EV['AUC','Value']

} #closes replicates


#### FIT FULL MODEL #####

SpData<-SpData[complete.cases(SpData),]

if(equalPrevalence==TRUE) {
tw<-table(SpData$PA)#n of pres and abs
wp<-(sum(tw)/2)/tw['1'] #weight pres
wb<-(sum(tw)/2)/tw['0'] #weight abs
weights<-ifelse(SpData$PA==1, wp, wb)
} else {weights=rep(1, nrow(SpData))}

if(Model == 'GLM') {
	#F<-as.formula(paste0('PA ~ ', paste(names(training[,-c(1:3)]), collapse='+')))#only linear
 	F<-as.formula(paste0('PA ~ ', paste(paste0('poly(',names(SpData[,-c(1:3)]),',2)'), collapse='+'))) #first order polynomial
	mod<-suppressWarnings(stepAIC(glm(F, family='binomial', weights=weights, data=SpData), trace=FALSE))
	prediction<-predict(mod, SpData, type='response')
	PresSuit<-predict(mod, EnvPres, type='response')
	FutSuit<-predict(mod, EnvFut, type='response')
}
if(Model == 'MaxEnt') {
	mod<-maxnet(p=SpData$PA, data=SpData[,-c(1:3)], maxnet.formula(p=SpData$PA, data=SpData[,-c(1:3)], classes="lq"))
	prediction<-predict(mod, SpData, clamp=TRUE, type=MXPredType)[,1]
	PresSuit<-predict(mod, EnvPres, clamp=TRUE, type=MXPredType)[,1]
	FutSuit<-predict(mod, EnvFut, clamp=TRUE, type=MXPredType)[,1]
}

if(Model == 'RF') {
	F<-as.formula(paste0('PA ~ ', paste(names(SpData[,-c(1:3)]), collapse='+')))
	SpDataRF<-SpData
	SpDataRF$PA<-as.factor(SpDataRF$PA)

    mod <- randomForest(formula = F,
                        data = SpDataRF,
                        ntree = 500,
                        mtry = round(sqrt(length(names(SpData[,-c(1:3)])))),
                        importance = FALSE,
                        norm.votes = TRUE,
                        strata = factor(c(0,1))
                        )

 	prediction<-predict(mod, SpData, type='prob')[,'1']
 	PresSuit<-predict(mod, EnvPres, type='prob')[,'1']
	FutSuit<-predict(mod, EnvFut, type='prob')[,'1']
}

#AVERAGE VALUES ESTIMATED FROM REPLICATES
best_thresTSS_cv<-mean(tss_tr_cv, na.rm=TRUE)
best_thresEqSS_cv<-mean(eqss_tr_cv, na.rm=TRUE)

sens_cv<-mean(sens_cv, na.rm=TRUE)
spec_cv<-mean(spec_cv, na.rm=TRUE)
TSS_cv<-mean(tss_cv, na.rm=TRUE)
AUC_cv<-mean(auc_cv, na.rm=TRUE)
AUC_cvBin<-mean(auc_cvBin, na.rm=TRUE)

#VALIDATON ON "REAL" DATA

	PAcrop<-crop(PA, Ext_Rls) #validate only in the original extent (the area of interest)
	PAreal<-as.data.frame(rasterToPoints(PAcrop))

  	#only retain cells that match (this step is necessary because in same cases one or two cells mismatch and the function stops)
  	tmpPR<-EnvPres; tmpPR$Suit<-PresSuit; tmpPR$xy<-paste0(tmpPR$x,'_', tmpPR$y); PAreal$xy<-paste0(PAreal$x,'_', PAreal$y); tmpPR<-merge(tmpPR[,c('xy','x','y','Suit')], PAreal[,c('xy','layer')], by='xy')
  	names(SUIT_mask_pp)<-c('x','y','layer2'); SUIT_mask_pp$xy<-paste0(SUIT_mask_pp$x,'_', SUIT_mask_pp$y); tmpPR<-merge(tmpPR[,c('xy','x','y','Suit','layer')], SUIT_mask_pp[,c('xy','layer2')], by='xy')

	PA_tss_pres<-ifelse(tmpPR$Suit > best_thresTSS_cv, 1, 0)
	PA_eqss_pres<-ifelse(tmpPR$Suit > best_thresEqSS_cv, 1, 0)

	EV<-modEval(predProb=tmpPR$Suit, predBin=PA_tss_pres, obs=tmpPR$layer)
	sens_pres<-EV['Sensitivity','Value']
	spec_pres<-EV['Specificity','Value']
	TSS_pres<-EV['TSS','Value']
	AUC_pres<-EV['AUC','Value']

	dumbCor_pres<-cor(tmpPR$layer2, tmpPR$Suit, method='spearman')

	PA2crop<-crop(PA2, Ext_Rls) #validate only in the original extent (the area of interest)	
	PA2real<-as.data.frame(rasterToPoints(PA2crop))

  	#only retain cells that match (this step is necessary because in same cases one or two cells mismatch and the function stops)
  	tmpFT<-EnvFut; tmpFT$Suit<-FutSuit; tmpFT$xy<-paste0(tmpFT$x,'_', tmpFT$y); PA2real$xy<-paste0(PA2real$x,'_', PA2real$y); tmpFT<-merge(tmpFT[,c('xy','x','y','Suit')], PA2real[,c('xy','layer')], by='xy')
  	names(SUIT_mask_pp2)<-c('x','y','layer2'); SUIT_mask_pp2$xy<-paste0(SUIT_mask_pp2$x,'_', SUIT_mask_pp2$y); tmpFT<-merge(tmpFT[,c('xy','x','y','Suit', 'layer')], SUIT_mask_pp2[,c('xy','layer2')], by='xy')

	PA_tss_fut<-ifelse(tmpFT$Suit > best_thresTSS_cv, 1, 0)
	PA_eqss_fut<-ifelse(tmpFT$Suit > best_thresEqSS_cv, 1, 0)

	EV<-modEval(predProb=tmpFT$Suit, predBin=PA_tss_fut, obs=tmpFT$layer)
	sens_fut<-EV['Sensitivity','Value']
	spec_fut<-EV['Specificity','Value']
	TSS_fut<-EV['TSS','Value']
	AUC_fut<-EV['AUC','Value']

	dumbCor_fut<-cor(tmpFT$layer2, tmpFT$Suit, method='spearman')

	#Accuracy of habitat loss and gain estimates (compare the real habitat loss and gain with the estimated ones using TSS and SOR)
	#real habitat loss
	PAstack<-stack(PA, PA2)
	PAstack<-crop(PAstack, Ext_Rls) #CROP BY THE AREA OF INTEREST
	PA_t<-as.data.frame(rasterToPoints(PAstack))
	names(PA_t)<-c('x','y','Pres','Fut')
	PA_t$HabitatLoss<-ifelse(PA_t$Pres>PA_t$Fut, 1, 0) #cells that are predicted to be lost using tss binarization
	PA_t$HabitatGain<-ifelse(PA_t$Pres<PA_t$Fut, 1, 0) #cells that are predicted to be lost using tss binarization
	PA_t$id<-paste0(PA_t$x, '_', PA_t$y)
	
	realHabLoss_per<-100*table(PA_t$HabitatLoss)['1']/table(PA_t$Pres)['1'] #potential loss in the full area (buffer included)
	realHabGain_per<-100*table(PA_t$HabitatGain)['1']/table(PA_t$Pres)['1'] #potential loss in the full area (buffer included)


	EstSuitPF<-merge(tmpPR[,c('xy','x','y','Suit')], tmpFT[,c('xy','Suit')], by='xy')
  	names(EstSuitPF)<-c('id','x','y','PresSuit','FutSuit')
  	EstSuitPF$ProbChangeLoss<-EstSuitPF$PresSuit-EstSuitPF$FutSuit; EstSuitPF$ProbChangeLoss<-ifelse(EstSuitPF$ProbChangeLoss<0, 0, EstSuitPF$ProbChangeLoss)
  	EstSuitPF$ProbChangeGain<-EstSuitPF$FutSuit-EstSuitPF$PresSuit; EstSuitPF$ProbChangeGain<-ifelse(EstSuitPF$ProbChangeGain<0, 0, EstSuitPF$ProbChangeGain)
	#tss
	EstSuitPF$Pres_tss<-ifelse(EstSuitPF$PresSuit > tss_cv[i], 1, 0)
	EstSuitPF$Fut_tss<-ifelse(EstSuitPF$FutSuit > tss_cv[i], 1, 0)
	EstSuitPF$HabitatLoss_tss<-ifelse(EstSuitPF$Pres_tss > EstSuitPF$Fut_tss, 1, 0) #cells that are predicted to be lost using tss binarization
	EstSuitPF$HabitatGain_tss<-ifelse(EstSuitPF$Pres_tss < EstSuitPF$Fut_tss, 1, 0) #cells that are predicted to be gained using tss binarization

	#comparison
	RealVsEst<-merge(PA_t, EstSuitPF[, c('id','HabitatLoss_tss', 'HabitatGain_tss', 'ProbChangeLoss','ProbChangeGain')], by='id')
	EV_HL<-modEval(predProb=RealVsEst$ProbChangeLoss, predBin=RealVsEst$HabitatLoss_tss, obs=RealVsEst$HabitatLoss)
	EV_HG<-modEval(predProb=RealVsEst$ProbChangeGain, predBin=RealVsEst$HabitatGain_tss, obs=RealVsEst$HabitatGain)
	sens_HabitatLoss<- EV_HL['Sensitivity','Value']
	spec_HabitatLoss<-EV_HL['Specificity','Value']
	sens_HabitatGain<-EV_HG['Sensitivity','Value']
	spec_HabitatGain<-EV_HG['Specificity','Value']
	TSS_HabitatLoss<-EV_HL['TSS','Value']
	TSS_HabitatGain<-EV_HG['TSS','Value']
	
	AUC_HabitatLoss<-EV_HL['AUC','Value']#independent of threshold
	AUC_HabitatGain<-EV_HG['AUC','Value']


EnvPres$Suit<-PresSuit
EnvFut$Suit<-FutSuit

#predict present
EnvPres$PAtss<-ifelse(EnvPres$Suit>best_thresTSS_cv, 1, 0)
EnvPres$PAauc<-ifelse(EnvPres$Suit>best_thresEqSS_cv, 1, 0)

predSuitpres<-rasterFromXYZ(EnvPres[,c('x','y','Suit')])
predPApres_tss<-rasterFromXYZ(EnvPres[,c('x','y','PAtss')])
predPApres_auc<-rasterFromXYZ(EnvPres[,c('x','y','PAauc')])

#predict future
EnvFut$PAtss<-ifelse(EnvFut$Suit>best_thresTSS_cv, 1, 0)
EnvFut$PAauc<-ifelse(EnvFut$Suit>best_thresEqSS_cv, 1, 0)

predSuitfut<-rasterFromXYZ(EnvFut[,c('x','y','Suit')])
predPAfut_tss<-rasterFromXYZ(EnvFut[,c('x','y','PAtss')])
predPAfut_auc<-rasterFromXYZ(EnvFut[,c('x','y','PAauc')])

#Estimate changes
#TSS
predPApres_tss<-crop(predPApres_tss, Ext_Rls) #CROP TO AREA OF INTEREST
predPAfut_tss<-crop(predPAfut_tss, Ext_Rls) #CROP TO AREA OF INTEREST

delta_est_tss<-predPApres_tss-predPAfut_tss
tmp<-as.data.frame(freq(delta_est_tss)); tmp<-tmp[complete.cases(tmp),]
loss_est_tss<-tmp[tmp$value==1,'count']; loss_est_tss<-ifelse(length(loss_est_tss)==0, 0, loss_est_tss)
gain_est_tss<-tmp[tmp$value==-1,'count']; gain_est_tss<-ifelse(length(gain_est_tss)==0, 0, gain_est_tss)
tmp<-as.data.frame(freq(predPApres_tss)); tmp<-tmp[complete.cases(tmp),]
AreaPresent_est_tss<-tmp[tmp$value==1,'count']; AreaPresent_est_tss<-ifelse(length(AreaPresent_est_tss)==0, 0, AreaPresent_est_tss)
tmp<-as.data.frame(freq(predPAfut_tss)); tmp<-tmp[complete.cases(tmp),]
AreaFuture_est_tss<-tmp[tmp$value==1,'count']; AreaFuture_est_tss<-ifelse(length(AreaFuture_est_tss)==0, 0, AreaFuture_est_tss)
gain_est_per_tss<-100*gain_est_tss/AreaPresent_est_tss
loss_est_per_tss<-100*loss_est_tss/AreaPresent_est_tss

#AUC
predPApres_auc<-crop(predPApres_auc, Ext_Rls) #CROP TO AREA OF INTEREST
predPAfut_auc<-crop(predPAfut_auc, Ext_Rls) #CROP TO AREA OF INTEREST

delta_est_auc<-predPApres_auc-predPAfut_auc
tmp<-as.data.frame(freq(delta_est_auc)); tmp<-tmp[complete.cases(tmp),]
loss_est_auc<-tmp[tmp$value==1,'count']; loss_est_auc<-ifelse(length(loss_est_auc)==0, 0, loss_est_auc)
gain_est_auc<-tmp[tmp$value==-1,'count']; gain_est_auc<-ifelse(length(gain_est_auc)==0, 0, gain_est_auc)
tmp<-as.data.frame(freq(predPApres_auc)); tmp<-tmp[complete.cases(tmp),]
AreaPresent_est_auc<-tmp[tmp$value==1,'count']; AreaPresent_est_auc<-ifelse(length(AreaPresent_est_auc)==0, 0, AreaPresent_est_auc)
tmp<-as.data.frame(freq(predPAfut_auc)); tmp<-tmp[complete.cases(tmp),]
AreaFuture_est_auc<-tmp[tmp$value==1,'count']; AreaFuture_est_auc<-ifelse(length(AreaFuture_est_auc)==0, 0, AreaFuture_est_auc)
gain_est_per_auc<-100*gain_est_auc/AreaPresent_est_auc
loss_est_per_auc<-100*loss_est_auc/AreaPresent_est_auc


stats=data.frame(Var=c('nPres', 'nBack', 'SpeciesPrevalence_stArea', 'SpeciesPrevalence_gExt', 'HyperVolume_pres', 'HyperVolume_abs', 'MESS_5', 'MESS_10', 'MESS_25', 'MESSland_0', 'MESSland_ave', 'propNonEq', 'sdBias', 'cvBias', 'realHabLoss_per', 'realHabGain_per', #Info Fitting
	'Sens_pres', 'Spec_pres', 'TSS_pres', 'Sens_fut', 'Spec_fut', 'TSS_fut', 'AUC_pres', 'AUC_fut', 'spCor_pres', 'spCor_fut', #info accuracy with "real" data #'TSS_tr', 'EqSS_tr', 
	'Sens_HabitatLoss', 'Spec_HabitatLoss', 'TSS_HabitatLoss', 'Sens_HabitatGain', 'Spec_HabitatGain', 'TSS_HabitatGain','AUC_HabitatLoss','AUC_HabitatGain',  #info accuracy with "real" data
	'Sens_cv', 'Spec_cv', 'TSS_cv', 'TSS_tr_cv', 'AUC_cv', 'EqSS_tr_cv', #info accuracy by cross-validation
	'AreaPresent_tss', 'AreaFuture_tss', 'gain_tss', 'gain_per_tss', 'loss_tss', 'loss_per_tss', #area estimates using TSS
	'AreaPresent_auc', 'AreaFuture_auc', 'gain_auc', 'gain_per_auc', 'loss_auc', 'loss_per_auc'), #area estimates using AUC
	Value=c(nrow(pres), nrow(background), SpPrev, pa_ratio, HyperVolume1, HyperVolume0, MESS_5, MESS_10, MESS_25, PresFutSimilarity_propMoreThan0, PresFutSimilarity_Mean, propNonEq, sdBias, cvBias, realHabLoss_per, realHabGain_per,
		sens_pres, spec_pres, TSS_pres, sens_fut, spec_fut, TSS_fut, AUC_pres, AUC_fut, dumbCor_pres, dumbCor_fut, #best_thresTSS, best_thresEqSS, 
		sens_HabitatLoss, spec_HabitatLoss, TSS_HabitatLoss, sens_HabitatGain, spec_HabitatGain, TSS_HabitatGain, AUC_HabitatLoss, AUC_HabitatGain,
		sens_cv, spec_cv, TSS_cv, best_thresTSS_cv, AUC_cv, best_thresEqSS_cv, 
		AreaPresent_est_tss, AreaFuture_est_tss, gain_est_tss, gain_est_per_tss, loss_est_tss, loss_est_per_tss,
		AreaPresent_est_auc, AreaFuture_est_auc, gain_est_auc, gain_est_per_auc, loss_est_auc, loss_est_per_auc)
	)

if(Plot==TRUE) {
par(mfrow=c(4,2))
plot(suit, main='Present suitability')
plot(PA, main='Present suitability')
plot(suit2, main='Future suitability')
plot(PA2, main='Future Presence')
plot(predSuitpres, main='Present suitability (estimated)')
plot(predPApres_tss, main='Present suitability (estimated)')
plot(predSuitfut, main='Future suitability (estimated)')
plot(predPAfut_tss, main='Future Presence (estimated)')
}

OUT<-list(predSuitpres=predSuitpres, predSuitfut=predSuitfut, 
	predPApres_tss=predPApres_tss, predPAfut_tss=predPAfut_tss,
	predPApres_auc=predPApres_auc, predPAfut_tss=predPAfut_auc,
	stats=stats, varDropped=excl)

return(OUT)

}


modEval<-function(predProb, predBin, obs){

	Rate<-ifelse(predBin==1 & obs==1, 'TruePositive', NA)
	Rate<-ifelse(predBin==0 & obs==0, 'TrueNegative', Rate)
	Rate<-ifelse(predBin==0 & obs==1, 'FalseNegative', Rate)
	Rate<-ifelse(predBin==1 & obs==0, 'FalsePositive', Rate)

	TP=table(Rate)['TruePositive']; TP<-ifelse(is.na(TP), 0, TP); names(TP)<-'TP'
	TN=table(Rate)['TrueNegative']; TN<-ifelse(is.na(TN), 0, TN); names(TN)<-'TN'
	FP=table(Rate)['FalsePositive']; FP<-ifelse(is.na(FP), 0, FP); names(FP)<-'FP'
	FN=table(Rate)['FalseNegative']; FN<-ifelse(is.na(FN), 0, FN); names(FN)<-'FN'

	Sensitivity=TP/(TP+FN); if(is.na(Sensitivity)){Sensitivity=0}
	Specificity=TN/(TN+FP); if(is.na(Specificity)){Specificity=0}
	TSS=Sensitivity + Specificity - 1
	if(!is.null(predProb) & length(unique(obs))!=1){
	AUC<-roc(obs, predProb, direction="<", levels = c(0,1), quiet=TRUE)$auc
	} else {AUC<-NA}

out<-data.frame(row.names=c('Sensitivity', 'Specificity', 'TSS', 'AUC'), Value=c(Sensitivity, Specificity, TSS, AUC))
return(out)
}

findOptimum<-function(predProb, obs, n_tr=100) {

precs<-seq(min(predProb), max(predProb), length.out=n_tr)

OUT=data.frame(Treshold=rep(NA, length(precs)), Sensitivity=NA, Specificity=NA, TSS=NA, SSAbsDif=NA)

for (ti in 1:length(precs)) {
t=precs[ti]
pred_T=ifelse(predProb<t,0, 1)

ev<-modEval(predBin=pred_T, predProb=predProb, obs=obs)

Sensitivity<-ev['Sensitivity','Value']
Specificity<-ev['Specificity','Value']
TSS<-ev['TSS','Value']
SSAbsDif<-abs(Sensitivity - Specificity) #to calc equal sens - spec

OUT[ti,]<-data.frame(t, Sensitivity, Specificity, TSS, SSAbsDif)
}
return(OUT)
}
