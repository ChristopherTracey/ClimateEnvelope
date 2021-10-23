##########################################################################################################################################
#                                                                                                                                        #
# Author: XXXX XXXX                                                                                                                	     #
# e-mail: XXXXXX@XXXXX                                                                                                    				 #
# Reference:                                                                                                                             #
# XXXX et al. 2019. Assessing the reliability of species distribution projections in climate change research. XXXXXX                     #
#                                                                                                                                        #
##########################################################################################################################################

#Plot partial response curves of random forest analysis

plot3per9Panels<-function(
	models=c('TSS_cv', 'TSS_pres', 'TSS_fut'), Title='9 panels', Main='', Yaxis='TSS', Ylegend=c('TSS cv', 'TSS (presence)', 'TSS (future)'),
	col1=rgb(27/255,158/255,119/255), colCI1=rgb(27/255,158/255,119/255, 0.5), col2=rgb(117/255,112/255,179/255), colCI2=rgb(117/255,112/255,179/255, 0.5), col3=rgb(217/255,95/255,2/255), colCI3=rgb(217/255,95/255,2/255,0.5), 
	linw=1.5, LoessSpan=0.5, distYlab=0.35, roundAxis=2, 
	dir='') {

#models=c('TSS_cv', 'TSS_pres', 'TSS_fut'); Title='9 panels'; Yaxis='TSS'; Ylegend=c('TSS cv', 'TSS (presence)', 'TSS (future)');
#col1=rgb(27/255,158/255,119/255); colCI1=rgb(27/255,158/255,119/255, 0.5); col2=rgb(117/255,112/255,179/255); colCI2=rgb(117/255,112/255,179/255, 0.5); col3=rgb(217/255,95/255,2/255); colCI3=rgb(217/255,95/255,2/255,0.5)
#linw=1.5; LoessSpan=0.5; distYlab=0.35; roundAxis=2
#dir='~/Documents/Progetti/VirtualSDMclimateProjection/SURFSARA_STUFF/figs/Combined/'

predictors<-c('SpeciesPrevalence_gExt', 'SamplePrevalence', 'nPresSampled', 'MESSland_ave', 'Bias', 'propNonEq','bufferP', 'TruePredictors', 'FalsePredictors')#'HyperVolume', 
Xlabels<-c('Species Prevalence', 'Sample Prevalence', 'n Presences', 'Env. similarity', 'Env. gradient sampled', 'Niche Filling','Geographic Extent', 'Relevant Predictors', 'Irrelevant Predictors')#'HyperVolume', 
#Ylabels<-c('real - est. range expansion (SOR)', 'real - est. range contraction (SOR)')


#MOD1
	mod1<-models[1]
	YHAT1<-get(paste0(mod1,'_YHAT'))
	RESP1<-get(paste0(mod1,'_RESP'))
	RESP_mn1<-apply(simplify2array(RESP1), 1:2, mean)
	
	LOESS_list1<-list()
	LOESS<-matrix(NA, nrow=nrow(RESP_mn1), ncol=ncol(RESP_mn1), dimnames=list(rownames(RESP_mn1), colnames(RESP_mn1)))
	for (k in 1:length(RESP1)) {LOESS_list1[[k]]<-LOESS}

	for (L in 1:length(RESP1)) { #for each species
	print(mod1)
	for (j in 1:length(predictors)){
		pred<-predictors[j]
		print(pred)
		l <-loess(RESP1[[L]][,pred] ~ YHAT1[,pred], span=LoessSpan)
		LOESS_list1[[L]][,pred]<-predict(l)
	}
	}
	RESP_mn1<-apply(simplify2array(LOESS_list1), 1:2, mean)
	RESP_se1<-apply(simplify2array(LOESS_list1), 1:2, sd)/sqrt(length(LOESS_list1))
	RESP_95CI1<-1.96*RESP_se1

#MOD 2
	mod2<-models[2]
	YHAT2<-get(paste0(mod2,'_YHAT'))
	RESP2<-get(paste0(mod2,'_RESP'))
	RESP_mn2<-apply(simplify2array(RESP2), 1:2, mean)

	LOESS_list2<-list()
	LOESS<-matrix(NA, nrow=nrow(RESP_mn2), ncol=ncol(RESP_mn2), dimnames=list(rownames(RESP_mn2), colnames(RESP_mn2)))
	for (k in 1:length(RESP2)) {LOESS_list2[[k]]<-LOESS}

	for (L in 1:length(RESP2)) { #for each species
	print(mod2)
	for (j in 1:length(predictors)){
		pred<-predictors[j]
		print(pred)
		l <-loess(RESP2[[L]][,pred] ~ YHAT2[,pred], span=LoessSpan)
		LOESS_list2[[L]][,pred]<-predict(l)
	}
	}
	RESP_mn2<-apply(simplify2array(LOESS_list2), 1:2, mean)
	RESP_se2<-apply(simplify2array(LOESS_list2), 1:2, sd)/sqrt(length(LOESS_list2))
	RESP_95CI2<-1.96*RESP_se2

#MOD 3
	mod3<-models[3]
	YHAT3<-get(paste0(mod3,'_YHAT'))
	RESP3<-get(paste0(mod3,'_RESP'))
	RESP_mn3<-apply(simplify2array(RESP3), 1:2, mean)

	LOESS_list3<-list()
	LOESS<-matrix(NA, nrow=nrow(RESP_mn3), ncol=ncol(RESP_mn3), dimnames=list(rownames(RESP_mn3), colnames(RESP_mn3)))
	for (k in 1:length(RESP3)) {LOESS_list3[[k]]<-LOESS}

	for (L in 1:length(RESP3)) { #for each species
	print(mod2)
	for (j in 1:length(predictors)){
		pred<-predictors[j]
		print(pred)
		l <-loess(RESP3[[L]][,pred] ~ YHAT3[,pred], span=LoessSpan)
		LOESS_list3[[L]][,pred]<-predict(l)
	}
	}
	RESP_mn3<-apply(simplify2array(LOESS_list3), 1:2, mean)
	RESP_se3<-apply(simplify2array(LOESS_list3), 1:2, sd)/sqrt(length(LOESS_list3))
	RESP_95CI3<-1.96*RESP_se3

	#range ylm
	ylm1<-range(c(RESP_mn1[,'SpeciesPrevalence_gExt'], RESP_mn1[,'SamplePrevalence'], RESP_mn1[,'nPresSampled'], RESP_mn1[,'MESSland_ave'], RESP_mn1[,'Bias'], RESP_mn1[,'bufferP'], RESP_mn1[,'propNonEq'], RESP_mn1[,'TruePredictors'], RESP_mn1[,'FalsePredictors']))#get(paste0(mod, '_', 'HyperVolume'))[,3], 
	ylm2<-range(c(RESP_mn2[,'SpeciesPrevalence_gExt'], RESP_mn2[,'SamplePrevalence'], RESP_mn2[,'nPresSampled'], RESP_mn2[,'MESSland_ave'], RESP_mn2[,'Bias'], RESP_mn2[,'bufferP'], RESP_mn2[,'propNonEq'], RESP_mn2[,'TruePredictors'], RESP_mn2[,'FalsePredictors']))#get(paste0(mod, '_', 'HyperVolume'))[,3], 
	ylm3<-range(c(RESP_mn3[,'SpeciesPrevalence_gExt'], RESP_mn3[,'SamplePrevalence'], RESP_mn3[,'nPresSampled'], RESP_mn3[,'MESSland_ave'], RESP_mn3[,'Bias'], RESP_mn3[,'bufferP'], RESP_mn3[,'propNonEq'], RESP_mn3[,'TruePredictors'], RESP_mn3[,'FalsePredictors']))#get(paste0(mod, '_', 'HyperVolume'))[,3], 
	#ylm3=c(0.5,0.7)
	ylm=range(c(ylm1, ylm2, ylm3))

	#save
#	pdf(paste0(dir, MODEL, '_', Title,'_PP_average.pdf'), width=5.5, height=16)
	pdf(paste0(dir, MODEL, '_', Title,'_PP_average.pdf'), width=16, height=6)
	
	#margins
	m1=3.5; m2=4; m3=3; m4=1

	par(mfrow=c(3,9), mar=c(m1, m2, m3, m4), las=1)	
#	par(mfrow=c(8,3), mar=c(m1, m2, m3, m4), las=1)

#CV
	for (i in 1:length(predictors)){
	
	pred<-predictors[i]
	par(mar=c(m1, m2, m3, m4), las=1)
	dif<-max(YHAT1[,pred]) - min(YHAT1[,pred])
	mar1<-min(YHAT1[,pred])-dif*distYlab; #print(mar1)
	mar2<-max(YHAT1[,pred])+dif*distYlab; #print(mar2)
	
	plot(YHAT1[,pred], RESP_mn1[,pred], type='l', xlab='', ylab=Yaxis, ylim=ylm1, col=col1, lwd=linw)
	polygon(x=c(YHAT1[,pred], rev(YHAT1[,pred])), y=c(RESP_mn1[,pred]-RESP_95CI1[,pred], rev(RESP_mn1[,pred]+RESP_95CI1[,pred])), col=colCI1, border=NA)
	mtext(text=Xlabels[i], side=1, line=2.2, cex=0.7)
		
	}
#PRES
	for (i in 1:length(predictors)){
	
	pred<-predictors[i]
	par(mar=c(m1, m2, m3, m4), las=1)
	dif<-max(YHAT1[,pred]) - min(YHAT1[,pred])
	mar1<-min(YHAT1[,pred])-dif*distYlab; #print(mar1)
	mar2<-max(YHAT1[,pred])+dif*distYlab; #print(mar2)
	
	plot(YHAT2[,pred], RESP_mn2[,pred], type='l', xlab='', ylab=Yaxis, ylim=ylm2, col=col2, lwd=linw)
	polygon(x=c(YHAT2[,pred], rev(YHAT2[,pred])), y=c(RESP_mn2[,pred]-RESP_95CI2[,pred], rev(RESP_mn2[,pred]+RESP_95CI2[,pred])), col=colCI2, border=NA)
	mtext(text=Xlabels[i], side=1, line=2.2, cex=0.7)
	
	}

#FUT
	for (i in 1:length(predictors)){
	
	pred<-predictors[i]
	par(mar=c(m1, m2, m3, m4), las=1)
	dif<-max(YHAT1[,pred]) - min(YHAT1[,pred])
	mar1<-min(YHAT1[,pred])-dif*distYlab; #print(mar1)
	mar2<-max(YHAT1[,pred])+dif*distYlab; #print(mar2)
	
	plot(YHAT3[,pred], RESP_mn3[,pred], type='l', xlab='', ylab=Yaxis, ylim=ylm3, col=col3, lwd=linw)
	polygon(x=c(YHAT3[,pred], rev(YHAT3[,pred])), y=c(RESP_mn3[,pred]-RESP_95CI3[,pred], rev(RESP_mn3[,pred]+RESP_95CI3[,pred])), col=colCI3, border=NA)
	mtext(text=Xlabels[i], side=1, line=2.2, cex=0.7)
		
	}
	mtext(text=Main, side=3, line=-2, outer=TRUE, font=2)

	dev.off()

	#plot legend
	pdf(paste0(dir, MODEL, '_', Title,'_legend.pdf'), width=8, height=4)
	plot(1,1, col='white', axes=FALSE, xlim=c(-2,4),ylim=c(-2,1), xlab='', ylab='')
	lines(x=c(-2, -1), y=c(0.5, 0.5), col=col1, lwd=5)
	lines(x=c(-2, -1), y=c(-0.5, -0.5), col=col2, lwd=5)
	lines(x=c(-2, -1), y=c(-1.5, -1.5), col=col3, lwd=5)
	text(x=-0.7, y=0.5, labels=Ylegend[1], cex=2.5, pos=4)
	text(x=-0.7, y=-0.5, labels=Ylegend[2], cex=2.5, pos=4)
	text(x=-0.7, y=-1.5, labels=Ylegend[3], cex=2.5, pos=4)
	dev.off()
}

MODEL='GLM'
load(paste0('partialResponses_',MODEL,'.Rdata'))

plot3per9Panels(models=c('TSS_cv', 'TSS_pres', 'TSS_fut'), Title='', Main='', Yaxis='TSS', Ylegend=c('TSS cv', 'TSS (presence)', 'TSS (future)'))

load(paste0('partialResponses_',MODEL,'_AUC.Rdata'))

plot3per9Panels(models=c('AUC_cv', 'AUC_pres', 'AUC_fut'), Title='', Main='', Yaxis='AUC', Ylegend=c('AUC cv', 'AUC (presence)', 'AUC (future)'))



plot2per9Panels<-function(
	models=c('TSS_cv', 'TSS_pres'), Title='2per9 panels', Main='', Yaxis1='TSS', Yaxis2='TSS', Ylegend=c('TSS cv', 'TSS (presence)'),
	col1=rgb(117/255,112/255,179/255), colCI1=rgb(117/255,112/255,179/255, 0.5), col2=rgb(217/255,95/255,2/255), colCI2=rgb(217/255,95/255,2/255,0.5), 
	linw=1.5, LoessSpan=0.5, distYlab=0.35, roundAxis=2, 
	dir='~/Documents/Progetti/VirtualSDMclimateProjection/Codes/R1/prove/res/') {

#models=c('TSS_cv', 'TSS_pres', 'TSS_fut'); Title='9 panels'; Yaxis='TSS'; Ylegend=c('TSS cv', 'TSS (presence)', 'TSS (future)');
#col1=rgb(27/255,158/255,119/255); colCI1=rgb(27/255,158/255,119/255, 0.5); col2=rgb(117/255,112/255,179/255); colCI2=rgb(117/255,112/255,179/255, 0.5); col3=rgb(217/255,95/255,2/255); colCI3=rgb(217/255,95/255,2/255,0.5)
#linw=1.5; LoessSpan=0.5; distYlab=0.35; roundAxis=2
#dir='~/Documents/Progetti/VirtualSDMclimateProjection/SURFSARA_STUFF/figs/Combined/'

predictors<-c('SpeciesPrevalence_gExt', 'SamplePrevalence', 'nPresSampled', 'MESSland_ave', 'Bias', 'propNonEq','bufferP', 'TruePredictors', 'FalsePredictors')#'HyperVolume', 
Xlabels<-c('Species Prevalence', 'Sample Prevalence', 'n Presences', 'Env. similarity','Env. gradient sampled', 'Niche Filling','Geographic Extent', 'Relevant Predictors', 'Irrelevant Predictors')#'HyperVolume', 
#Ylabels<-c('real - est. range expansion (SOR)', 'real - est. range contraction (SOR)')


#MOD1
	mod1<-models[1]
	YHAT1<-get(paste0(mod1,'_YHAT'))
	RESP1<-get(paste0(mod1,'_RESP'))
	RESP_mn1<-apply(simplify2array(RESP1), 1:2, mean)
	
	LOESS_list1<-list()
	LOESS<-matrix(NA, nrow=nrow(RESP_mn1), ncol=ncol(RESP_mn1), dimnames=list(rownames(RESP_mn1), colnames(RESP_mn1)))
	for (k in 1:length(RESP1)) {LOESS_list1[[k]]<-LOESS}

	for (L in 1:length(RESP1)) { #for each species
	print(mod1)
	for (j in 1:length(predictors)){
		pred<-predictors[j]
		print(pred)
		l <-loess(RESP1[[L]][,pred] ~ YHAT1[,pred], span=LoessSpan)
		LOESS_list1[[L]][,pred]<-predict(l)
	}
	}
	RESP_mn1<-apply(simplify2array(LOESS_list1), 1:2, mean)
	RESP_se1<-apply(simplify2array(LOESS_list1), 1:2, sd)/sqrt(length(LOESS_list1))
	RESP_95CI1<-1.96*RESP_se1

#MOD 2
	mod2<-models[2]
	YHAT2<-get(paste0(mod2,'_YHAT'))
	RESP2<-get(paste0(mod2,'_RESP'))
	RESP_mn2<-apply(simplify2array(RESP2), 1:2, mean)

	LOESS_list2<-list()
	LOESS<-matrix(NA, nrow=nrow(RESP_mn2), ncol=ncol(RESP_mn2), dimnames=list(rownames(RESP_mn2), colnames(RESP_mn2)))
	for (k in 1:length(RESP2)) {LOESS_list2[[k]]<-LOESS}

	for (L in 1:length(RESP2)) { #for each species
	print(mod2)
	for (j in 1:length(predictors)){
		pred<-predictors[j]
		print(pred)
		l <-loess(RESP2[[L]][,pred] ~ YHAT2[,pred], span=LoessSpan)
		LOESS_list2[[L]][,pred]<-predict(l)
	}
	}
	RESP_mn2<-apply(simplify2array(LOESS_list2), 1:2, mean)
	RESP_se2<-apply(simplify2array(LOESS_list2), 1:2, sd)/sqrt(length(LOESS_list2))
	RESP_95CI2<-1.96*RESP_se2

	#range ylm
	ylm1<-range(c(RESP_mn1[,'SpeciesPrevalence_gExt'], RESP_mn1[,'SamplePrevalence'], RESP_mn1[,'nPresSampled'], RESP_mn1[,'MESSland_ave'], RESP_mn1[,'Bias'], RESP_mn1[,'bufferP'], RESP_mn1[,'propNonEq'], RESP_mn1[,'TruePredictors'], RESP_mn1[,'FalsePredictors']))#get(paste0(mod, '_', 'HyperVolume'))[,3], 
	ylm2<-range(c(RESP_mn2[,'SpeciesPrevalence_gExt'], RESP_mn2[,'SamplePrevalence'], RESP_mn2[,'nPresSampled'], RESP_mn2[,'MESSland_ave'], RESP_mn2[,'Bias'], RESP_mn2[,'bufferP'], RESP_mn2[,'propNonEq'], RESP_mn2[,'TruePredictors'], RESP_mn2[,'FalsePredictors']))#get(paste0(mod, '_', 'HyperVolume'))[,3], 

	#save
	pdf(paste0(dir, MODEL, '_', Title,'_PP_average.pdf'), width=16, height=4)
	
	#margins
	m1=4; m2=4; m3=3; m4=1
	
	par(mfrow=c(2,9), mar=c(m1, m2, m3, m4), las=1)
	for (i in 1:length(predictors)){
	
	pred<-predictors[i]
	par(mar=c(m1, m2, m3, m4), las=1)
	dif<-max(YHAT1[,pred]) - min(YHAT1[,pred])
	mar1<-min(YHAT1[,pred])-dif*distYlab; #print(mar1)
	mar2<-max(YHAT1[,pred])+dif*distYlab; #print(mar2)
	
	plot(YHAT1[,pred], RESP_mn1[,pred], type='l', xlab='', ylab=Yaxis1, ylim=ylm1, col=col1, lwd=linw)
	polygon(x=c(YHAT1[,pred], rev(YHAT1[,pred])), y=c(RESP_mn1[,pred]-RESP_95CI1[,pred], rev(RESP_mn1[,pred]+RESP_95CI1[,pred])), col=colCI1, border=NA)
	mtext(text=Xlabels[i], side=1, line=2.2, cex=0.7)
	
	}

	for (i in 1:length(predictors)){
	
	pred<-predictors[i]
	par(mar=c(m1, m2, m3, m4), las=1)
	dif<-max(YHAT1[,pred]) - min(YHAT1[,pred])
	mar1<-min(YHAT1[,pred])-dif*distYlab; #print(mar1)
	mar2<-max(YHAT1[,pred])+dif*distYlab; #print(mar2)
	
	plot(YHAT2[,pred], RESP_mn2[,pred], type='l', xlab='', ylab=Yaxis2, ylim=ylm2, col=col2, lwd=linw)
	polygon(x=c(YHAT2[,pred], rev(YHAT2[,pred])), y=c(RESP_mn2[,pred]-RESP_95CI2[,pred], rev(RESP_mn2[,pred]+RESP_95CI2[,pred])), col=colCI2, border=NA)
	mtext(text=Xlabels[i], side=1, line=2.2, cex=0.7)
	
	}
	mtext(text=Main, side=3, line=-2, outer=TRUE, font=2)

	dev.off()

	#plot legend
	pdf(paste0(dir, MODEL, '_', Title,'_legend.pdf'), width=8, height=4)
	plot(1,1, col='white', axes=FALSE, xlim=c(-2,4),ylim=c(-1,1), xlab='', ylab='')
	lines(x=c(-2, -1), y=c(0.5, 0.5), col=col1, lwd=5)
	lines(x=c(-2, -1), y=c(-0.5, -0.5), col=col2, lwd=5)
	text(x=-0.7, y=0.5, labels=Ylegend[1], cex=2.5, pos=4)
	text(x=-0.7, y=-0.5, labels=Ylegend[2], cex=2.5, pos=4)
	dev.off()
}

MODEL='GLM'
load(paste0('/partialResponses_',MODEL,'.Rdata'))
plot2per9Panels(models=c('est_gain_per_tss', 'est_loss_per_tss'), Title='9 panels est gain vs est loss', Main='', Yaxis1='Range Expansion %', Yaxis2='Range Contraction %', Ylegend=c('Range Expansion %', 'Range Contraction %'))




