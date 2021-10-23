##########################################################################################################################################
#                                                                                                                                        #
# Author: XXXX XXXX                                                                                                                	     #
# e-mail: XXXXXX@XXXXX                                                                                                    				 #
# Reference:                                                                                                                             #
# XXXX et al. 2019. Assessing the reliability of species distribution projections in climate change research. XXXXXX                     #
#                                                                                                                                        #
##########################################################################################################################################

#This code uses random forest to estimate the effect and importance of different variables on the estimated and true TSS

library(ranger)
library(pdp)

MODEL='GLM'
data<-read.csv(paste0('OUTPUT_Simulation_',MODEL,'_R1.csv'))
data<-data[!is.na(data$SpID),]

data$delta_loss<-data$loss_per2-data$est_loss_per_tss
data$delta_gain<-data$gain_per2-data$est_gain_per_tss

res<-20
predictors<-c('SpeciesPrevalence_gExt', 'SamplePrevalence', 'nPresSampled', 'Bias', 'propNonEq','bufferP', 'TruePredictors', 'FalsePredictors','MESSland_ave')#'HyperVolume', 'MESS_5', 
responses<-c('Sens_pres', 'Spec_pres', 'TSS_pres', 'Sens_fut', 'Spec_fut', 'TSS_fut', 'TSS_HabitatLoss', 'TSS_HabitatGain', 'delta_gain', 'delta_loss', 'est_gain_per_tss', 'est_loss_per_tss', 'TSS_cv')

data<-data[complete.cases(data[,c(predictors, responses)]),]
data<-data[data$SamplePrevalence<=1,]

spids<-unique(data$SpID)

out_imp<-matrix(NA, ncol=length(predictors), nrow=length(spids))

#ALL
for (i_MOD in 1:length(responses)){

Y<-responses[i_MOD]
F<-formula(paste0(Y, ' ~ SpeciesPrevalence_gExt + SamplePrevalence + nPresSampled + propNonEq + Bias + bufferP + TruePredictors + FalsePredictors + MESSland_ave'))

	for (i in 1:length(spids)){
	
	spid<-spids[i]
	print(spid)
	Dsp=data[data$SpID==spid,]

	mod<-ranger(F, 
		data=Dsp, mtry=3, importance='permutation', write.forest=TRUE)
	
	#Importance
	imp<-importance(mod)
	out_imp[i,]<-100*imp/sum(imp)#rescale to 100

	#create matrices to fill
	yhat<-matrix(nrow=res, ncol=length(predictors), dimnames=list(1:res, predictors))
	resp<-matrix(nrow=res, ncol=length(predictors), dimnames=list(1:res, predictors))
	#loop over predictors and store partial responses
		for (j in 1:length(predictors)){
		print(paste0('pred ',j))
		pred<-predictors[j]
		pgrid<-data.frame(seq(min(data[,pred], na.rm=TRUE), max(data[,pred], na.rm=TRUE), length.out=res)); names(pgrid)<-pred
		p<-partial(mod, pred.var=pred, pred.grid=pgrid)
		yhat[,j]<-p[,1]
		resp[,j]<-p[,2]
	}
	assign(paste0('yhat',i), yhat)
	assign(paste0('resp',i), resp)
	}

	#average Importances
	ave_imp<-apply(out_imp, 2, mean, na.rm=TRUE)
	names(ave_imp)<-names(importance(mod))
	assign(paste0(Y,'_Imp_mn'), ave_imp)

	sd_imp<-apply(out_imp, 2, sd, na.rm=TRUE)
	names(sd_imp)<-names(importance(mod))
	assign(paste0(Y,'_Imp_sd'), sd_imp)

	#average partial responses over species
	YHAT<-list()
	RESP<-list()
	for (i in 1:length(spids)){
		YHAT[[i]]<-get(paste0('yhat',i))
		RESP[[i]]<-get(paste0('resp',i))
	}
	YHAT_mn<-apply(simplify2array(YHAT), 1:2, mean)
	RESP_mn<-apply(simplify2array(RESP), 1:2, mean)
	RESP_sd<-apply(simplify2array(RESP), 1:2, sd)

	assign(paste0(Y, '_YHAT'), YHAT_mn)
	assign(paste0(Y, '_RESP'), RESP)
	assign(paste0(Y, '_RESP_mn'), RESP_mn)
	assign(paste0(Y, '_RESP_sd'), RESP_sd)
}

save(Sens_pres_RESP, Spec_pres_RESP, TSS_pres_RESP, Sens_fut_RESP, Spec_fut_RESP, TSS_fut_RESP, TSS_HabitatLoss_RESP, TSS_HabitatGain_RESP, delta_gain_RESP, delta_loss_RESP, est_gain_per_tss_RESP, est_loss_per_tss_RESP, TSS_cv_RESP,
	Sens_pres_RESP_mn, Spec_pres_RESP_mn, TSS_pres_RESP_mn, Sens_fut_RESP_mn, Spec_fut_RESP_mn, TSS_fut_RESP_mn, TSS_HabitatLoss_RESP_mn, TSS_HabitatGain_RESP_mn, delta_gain_RESP_mn, delta_loss_RESP_mn, est_gain_per_tss_RESP_mn, est_loss_per_tss_RESP_mn, TSS_cv_RESP_mn, 
	Sens_pres_RESP_sd, Spec_pres_RESP_sd, TSS_pres_RESP_sd, Sens_fut_RESP_sd, Spec_fut_RESP_sd, TSS_fut_RESP_sd, TSS_HabitatLoss_RESP_sd, TSS_HabitatGain_RESP_sd, delta_gain_RESP_sd, delta_loss_RESP_sd, est_gain_per_tss_RESP_sd, est_loss_per_tss_RESP_sd, TSS_cv_RESP_sd, 
	Sens_pres_YHAT, Spec_pres_YHAT, TSS_pres_YHAT, Sens_fut_YHAT, Spec_fut_YHAT, TSS_fut_YHAT, TSS_HabitatLoss_YHAT, TSS_HabitatGain_YHAT, delta_gain_YHAT, delta_loss_YHAT, est_gain_per_tss_YHAT, est_loss_per_tss_YHAT, TSS_cv_YHAT, 
	file=paste0('partialResponses_',MODEL,'_R1.Rdata'))

save(Sens_pres_Imp_mn, Spec_pres_Imp_mn, TSS_pres_Imp_mn, Sens_fut_Imp_mn, Spec_fut_Imp_mn, TSS_fut_Imp_mn, TSS_HabitatLoss_Imp_mn, TSS_HabitatGain_Imp_mn, delta_gain_Imp_mn, delta_loss_Imp_mn, est_gain_per_tss_Imp_mn, est_loss_per_tss_Imp_mn, TSS_cv_Imp_mn, 
	Sens_pres_Imp_sd, Spec_pres_Imp_sd, TSS_pres_Imp_sd, Sens_fut_Imp_sd, Spec_fut_Imp_sd, TSS_fut_Imp_sd, TSS_HabitatLoss_Imp_sd, TSS_HabitatGain_Imp_sd, delta_gain_Imp_sd, delta_loss_Imp_sd, est_gain_per_tss_Imp_sd, est_loss_per_tss_Imp_sd, TSS_cv_Imp_sd, 
	file=paste0('Importances_',MODEL,'_R1.Rdata'))




