##########################################################################################################################################
#                                                                                                                                        #
# Author: XXXX XXXX                                                                                                                	     #
# e-mail: XXXXXX@XXXXX                                                                                                    				 #
# Reference:                                                                                                                             #
# XXXX et al. 2019. Assessing the reliability of species distribution projections in climate change research. XXXXXX                     #
#                                                                                                                                        #
##########################################################################################################################################

#This code generates virtual species and fit species distribution models (SDMs) using diverse settings
#It then validate the SDMs using a split sample validation, and against the virtual reality for both present and future predictions
#the code is structured as an array job

library(virtualspecies)
library(PresenceAbsence)
library(maxnet)
library(dismo)
library(rgeos)
library(hypervolume)
library(clhs)
library(raster)
library(rgdal)

source('Functions.R', chdir = TRUE)

#set grid to resample rasters
res=0.1
mat<-matrix(0, ncol=360/res, nrow=180/res)
grid<-raster(mat, xmn=-180, xmx=180, ymn=-90, ymx=90, crs=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

#load rasters
dirPres<-('bioclim_Chelsa_current/')
dirFut<-('Bioclim_Chelsa_2050_RCP85/')

filesPres<-dir(dirPres, full.names=TRUE)
filesFut<-dir(dirFut, full.names=TRUE)

s<-stack(filesPres)
s2<-stack(filesFut)

#load raster with land blocks
LandBlocks <- raster('/vol/milkunB/lsantini/Vspecies/envVar/LandBlocks.tif')
LandBlocks<-resample(LandBlocks, grid, method='ngb')
EU=87
AF=2
AU=4
NAm=1
SA=3

#load human influence index map
HIImap<-raster('/vol/milkunB/lsantini/Vspecies/envVar/HFP2009_ll.tif')
HIImap<-resample(HIImap, grid)

#set treatments
presSeq<-c(10, 25, 50, 100, 250, 500, 1000)
buffP<-c(0, 100, 500, 5000, 50000)
BiasQ<-c(33, 66, 100)
nonEqT<-c(0.33, 0.66, 1)
sPrevSeq<-c(0.01, 0.1, 1)
nTruePredictors<-c(0, 3, 6) #now is always 4 variables used. What Happens if the number of true stays the same but false increase?
nFalsePredictors<-c(0, 3, 6) #now is always 4 variables used. What Happens if the number of true stays the same but false increase?

PresAbs<-expand.grid(presSeq, buffP, BiasQ, nonEqT, sPrevSeq, nTruePredictors, nFalsePredictors) #combinations
names(PresAbs)<-c('nPres','buffP', 'Bias', 'nonEqT', 'sPrev', 'nTruePredictors', 'nFalsePredictors')
PresAbs<-PresAbs[!(PresAbs$nTruePredictors==0 & PresAbs$nFalsePredictors==0),] #exclude cases with no variables
PresAbs$nAbs<-PresAbs$nPres/PresAbs$sPrev

#sample from multivariate distribution using a conditioned latin hypercube to maximize the range of variables
set.seed(1)
index <- clhs(PresAbs, size = 500, progress = FALSE, iter = 1000) 
PresAbs<-PresAbs[index,]

MODEL='GLM'

nrows=nrow(PresAbs)#*nSp
OUT<-data.frame(SpID=NA, AreaPresent=NA, AreaFuture=NA, loss_per=NA, gain_per=NA, loss_per2=NA, gain_per2=NA, HyperVolume=NA, SpeciesPrevalence_stArea=NA, SpeciesPrevalence_gExt=NA, Continent=NA, PropRealized=NA,
				nPres=rep(NA, nrows), nAbs=rep(NA, nrows), bufferP=rep(NA, nrows), 
				nPresSampled=NA, nAbsSampled=NA, SamplePrevalence=NA, propNonEq=NA, Bias=NA, sdBias=NA, cvBias=NA, MESS_5=NA, MESS_10=NA, MESS_25=NA,
				HyperVolume_presence=NA, HyperVolume_absence=NA, Sens_pres=NA, Spec_pres=NA, TSS_pres=NA, Sens_fut=NA, Spec_fut=NA, TSS_fut=NA, AUC_pres=NA, AUC_fut=NA, #TSS_tr=NA, EqSS_tr=NA, 
				SensHabitatLoss=NA, SpecHabitatLoss=NA, TSS_HabitatLoss=NA, SensHabitatGain=NA, SpecHabitatGain=NA, TSS_HabitatGain=NA, AUC_HabitatLoss=NA, AUC_HabitatGain=NA,
				EstAreaPresent_tss=NA, EstAreaFuture_tss=NA, est_loss_per_tss=NA, est_gain_per_tss=NA, Sens_cv=NA, Spec_cv=NA, TSS_cv=NA, TSS_tr_cv=NA, 
				EstAreaPresent_auc=NA, EstAreaFuture_auc=NA, est_loss_per_auc=NA, est_gain_per_auc=NA, AUC_cv=NA, EqSS_tr_cv=NA 
				)
row=0

i <- as.numeric(commandArgs(trailingOnly = TRUE))

set.seed(i)

	#generate new random species
	repeat{
		VS<-try(GenerateRandomSP_19bioclim(LandBlocks=LandBlocks, nImpVar=6, RandomExtent=TRUE, SamplingExtent=5, stack_pres=s, stack_future=s2, stack_PCA=sPCA, n=100, minCell=5, randomContinent=TRUE, CONT=87, HyperVolumeNiche=TRUE, HypervolumeGaussian=FALSE, PAlogistic=FALSE, Plot=FALSE))
		if(class(VS)!='try-error') {break}
		}


	stats1<-VS[['stats']]
	ST<-VS[['suit_pt']]
	ST2<-VS[['suit2_pt']]
	PA<-VS[['PA_pt']]
	PA2<-VS[['PA2_pt']]
	Ext_Rls<-VS[['PA_rl']] #PA crop to the realized range (used for sampling, not for testing using real data)
 	cellsRealiz<-table(rasterToPoints(Ext_Rls)[,3])['1'] #n cells of realized occupancy
 	cellsPotential<-table(rasterToPoints(PA)[,3])['1'] #n cells of potential occupancy in the continent (it tell us how much information on the niche is missed that does not depend on the modeling setting)
 	PropRealized<-cellsRealiz/cellsPotential #in reality this depends on the size of the continent...
 	LB<-VS[['LB']]
 	ImpVar<-VS[['ImpVar']]

	AreaPresent<-stats1[stats1$Var=='AreaPresent','Value']
	AreaFuture<-stats1[stats1$Var=='AreaFuture','Value']
	gain_per<-stats1[stats1$Var=='gain_per','Value']
	loss_per<-stats1[stats1$Var=='loss_per','Value']
	HyperVolume<-stats1[stats1$Var=='HyperVolume','Value']
	Continent<-as.numeric(na.omit(unique(values(LB))))


for (j in 1:nrow(PresAbs)) {

	row=row+1
	nPres=PresAbs[j,'nPres']
	nAbs=PresAbs[j,'nAbs']
	buffP=PresAbs[j,'buffP']
	BiasQ=PresAbs[j,'Bias']
	nonEqT=PresAbs[j,'nonEqT']
	nTruePredictors=PresAbs[j,'nTruePredictors']
	nFalsePredictors=PresAbs[j,'nFalsePredictors']

print(paste(j, 'out of', nrow(PresAbs)))
print(paste('nPres =', nPres, '; nAbs =', nAbs, '; buffP =', buffP, '; Bias =', BiasQ, '; nonEq =', nonEqT, '; nTruePredictors =', nTruePredictors, '; nFalsePredictors =', nFalsePredictors))

	#select variables
	allVar<-c('bio_01', 'bio_02', 'bio_03', 'bio_04', 'bio_05', 'bio_06', 'bio_07', 'bio_08', 'bio_09', 'bio_10', 'bio_11', 'bio_12', 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18', 'bio_19')
	nonImpVar<-allVar[!allVar %in% ImpVar]
	varRet<-sample(ImpVar, nTruePredictors) #important variables used for fitting
	varAdd<-nFalsePredictors
	varAdd<-sample(nonImpVar, varAdd) #variables used for fitting but not important
	varFit<-c(varRet, varAdd) #variables used for fitting
	sb<-s[[varFit]]
	s2b<-s2[[varFit]]

	#Sample points and fit SDM
	SDM<-try(FitSDM(LB=LB, SUIT=ST, SUIT2=ST2, PA=PA, PA2=PA2, Ext_Rls=Ext_Rls, stack_pres=sb, stack_future=s2b, truePred=varRet, falsePred=varAdd, Model=MODEL, MXPredType='cloglog', ExcludeCollinear=TRUE, equalPrevalence=TRUE, n_pres=nPres, n_back=nAbs, backgroundWithinBuffer=TRUE, bufferP=buffP, SpatVal=TRUE, replicates=10, rangeSDM=FALSE, AbsNotInPres=FALSE, HyperVolumePresAbs=FALSE, HypervolumeGaussian=FALSE, Bias=TRUE, BiasQ=BiasQ, BiasMap=NULL, nonEquilibrium=TRUE, HIImap=HIImap, nonEquilibriumQuantileValue=nonEqT, Plot=FALSE))

	if(class(SDM)=='try-error') {next}

	stats2<-SDM[['stats']]

	varDropped<-SDM[['varDropped']] #which variables were dropped because collinear?
	varRet<-varRet[!varRet %in% varDropped] #remove those that were dropped
	varAdd<-varAdd[!varAdd %in% varDropped] #remove those that were dropped

#info species
	OUT[row, 'SpID']<-paste0('Sp',i)
	OUT[row, 'AreaPresent']<-AreaPresent
	OUT[row, 'AreaFuture']<-AreaFuture
	OUT[row, 'gain_per']<-gain_per #percentage considering the original study area
	OUT[row, 'loss_per']<-loss_per #percentage considering the original study area
	OUT[row, 'gain_per2']<-stats2[stats2$Var=='realHabGain_per','Value'] #percentage considering the potential presence and expansion/contraction in the expanded study area
	OUT[row, 'loss_per2']<-stats2[stats2$Var=='realHabLoss_per','Value'] #percentage considering the potential presence and expansion/contraction in the expanded study area
	OUT[row, 'HyperVolume']<-HyperVolume
	OUT[row, 'SpeciesPrevalence_stArea']<-stats2[stats2$Var=='SpeciesPrevalence_stArea','Value'] #not the original Species prevalence calculated in GenerateRanddomSp(), but that dependent on the study area calculated by FitSDM()
	OUT[row, 'SpeciesPrevalence_gExt']<-stats2[stats2$Var=='SpeciesPrevalence_gExt','Value'] #not the original Species prevalence calculated in GenerateRanddomSp(), but that dependent on the study area calculated by FitSDM()
	OUT[row, 'Continent']<-Continent
	OUT[row, 'PropRealized']<-PropRealized

#info settings
	OUT[row, 'nPres']<-nPres
	OUT[row, 'nAbs']<-nAbs
	OUT[row, 'bufferP']<-buffP
	OUT[row, 'nPresSampled']<-stats2[stats2$Var=='nPres','Value']
	OUT[row, 'nAbsSampled']<-stats2[stats2$Var=='nBack','Value']
	OUT[row, 'SamplePrevalence']<-OUT[row, 'nPresSampled'] / OUT[row, 'nAbsSampled']
	OUT[row, 'propNonEq']<-stats2[stats2$Var=='propNonEq','Value']
	OUT[row, 'Bias']<-BiasQ
	OUT[row, 'sdBias']<-stats2[stats2$Var=='sdBias','Value']
	OUT[row, 'cvBias']<-stats2[stats2$Var=='cvBias','Value']
	OUT[row, 'MESS_5']<-stats2[stats2$Var=='MESS_5','Value']
	OUT[row, 'MESS_10']<-stats2[stats2$Var=='MESS_10','Value']
	OUT[row, 'MESS_25']<-stats2[stats2$Var=='MESS_25','Value']
	OUT[row, 'MESSland_0']<-stats2[stats2$Var=='MESSland_0','Value']
	OUT[row, 'MESSland_ave']<-stats2[stats2$Var=='MESSland_ave','Value']
	OUT[row, 'TruePredictors']<-length(varRet)
	OUT[row, 'FalsePredictors']<-length(varAdd)

#info model
	OUT[row, 'HyperVolume_presence']<-stats2[stats2$Var=='HyperVolume_pres','Value']
	OUT[row, 'HyperVolume_absence']<-stats2[stats2$Var=='HyperVolume_abs','Value']
	OUT[row, 'Sens_pres']<-stats2[stats2$Var=='Sens_pres','Value']
	OUT[row, 'Spec_pres']<-stats2[stats2$Var=='Spec_pres','Value']
	OUT[row, 'TSS_pres']<-stats2[stats2$Var=='TSS_pres','Value']
	OUT[row, 'Sens_fut']<-stats2[stats2$Var=='Sens_fut','Value']
	OUT[row, 'Spec_fut']<-stats2[stats2$Var=='Spec_fut','Value']
	OUT[row, 'TSS_fut']<-stats2[stats2$Var=='TSS_fut','Value']
	OUT[row, 'AUC_pres']<-stats2[stats2$Var=='AUC_pres','Value']
	OUT[row, 'AUC_fut']<-stats2[stats2$Var=='AUC_fut','Value']
	OUT[row, 'spCor_pres']<-stats2[stats2$Var=='spCor_pres','Value']
	OUT[row, 'spCor_fut']<-stats2[stats2$Var=='spCor_fut','Value']

	OUT[row, 'Sens_HabitatLoss']<-stats2[stats2$Var=='Sens_HabitatLoss','Value']
	OUT[row, 'Spec_HabitatLoss']<-stats2[stats2$Var=='Spec_HabitatLoss','Value']
	OUT[row, 'TSS_HabitatLoss']<-stats2[stats2$Var=='TSS_HabitatLoss','Value']
	OUT[row, 'AUC_HabitatLoss']<-stats2[stats2$Var=='AUC_HabitatLoss','Value']

	OUT[row, 'Sens_HabitatGain']<-stats2[stats2$Var=='Sens_HabitatGain','Value']
	OUT[row, 'Spec_HabitatGain']<-stats2[stats2$Var=='Spec_HabitatGain','Value']
	OUT[row, 'TSS_HabitatGain']<-stats2[stats2$Var=='TSS_HabitatGain','Value']
	OUT[row, 'AUC_HabitatGain']<-stats2[stats2$Var=='AUC_HabitatGain','Value']

	OUT[row, 'EstAreaPresent_tss']<-stats2[stats2$Var=='AreaPresent_tss','Value']
	OUT[row, 'EstAreaFuture_tss']<-stats2[stats2$Var=='AreaFuture_tss','Value']
	OUT[row, 'est_gain_per_tss']<-stats2[stats2$Var=='gain_per_tss','Value']
	OUT[row, 'est_loss_per_tss']<-stats2[stats2$Var=='loss_per_tss','Value']
	OUT[row, 'Sens_cv']<-stats2[stats2$Var=='Sens_cv','Value']
	OUT[row, 'Spec_cv']<-stats2[stats2$Var=='Spec_cv','Value']
	OUT[row, 'TSS_cv']<-stats2[stats2$Var=='TSS_cv','Value']
	OUT[row, 'TSS_tr_cv']<-stats2[stats2$Var=='TSS_tr_cv','Value']
	
	OUT[row, 'EstAreaPresent_auc']<-stats2[stats2$Var=='AreaPresent_auc','Value']
	OUT[row, 'EstAreaFuture_auc']<-stats2[stats2$Var=='AreaFuture_auc','Value']
	OUT[row, 'est_gain_per_auc']<-stats2[stats2$Var=='gain_per_auc','Value']
	OUT[row, 'est_loss_per_auc']<-stats2[stats2$Var=='loss_per_auc','Value']
	OUT[row, 'AUC_cv']<-stats2[stats2$Var=='AUC_cv','Value']
	OUT[row, 'EqSS_tr_cv']<-stats2[stats2$Var=='EqSS_tr_cv','Value']

	}
write.csv(OUT, paste0('/OutArrayJob/OUT_',MODEL,'_Sp',i,'_nVariables_tr.csv'), row.names=FALSE)
