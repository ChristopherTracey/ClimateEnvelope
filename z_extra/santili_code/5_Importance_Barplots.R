##########################################################################################################################################
#                                                                                                                                        #
# Author: XXXX XXXX                                                                                                                      #
# e-mail: XXXXXX@XXXXX                                                                                                             #
# Reference:                                                                                                                             #
# XXXX et al. 2019. Assessing the reliability of species distribution projections in climate change research. XXXXXX                     #
#                                                                                                                                        #
##########################################################################################################################################

#plot importance estimates of random forest

library(ggplot2)
library(ggpubr)

ord<-rev(c('SpeciesPrevalence_gExt','SamplePrevalence','nPresSampled', 'MESSland_ave','Bias', 'propNonEq','bufferP','TruePredictors','FalsePredictors'))
labs<-rev(c('Species Prevalence', 'Sample Prevalence', 'n Presences', 'Env. similarity', 'Env. gradient sampled', 'Niche Filling', 'Geographic Extent', 'Relevant Predictors', 'Irrelevant Predictors'))

models<-c('RF', 'MaxEnt', 'GLM')
for (i in 1:length(models)) {
MODEL=models[i]
nSpecies=50
load(paste0('Importances_',MODEL,'_R1.Rdata'))

#TSS cv
imp_cv_tss<-data.frame(Var=names(TSS_cv_Imp_mn), Imp=as.numeric(TSS_cv_Imp_mn), se=as.numeric(TSS_cv_Imp_sd)/sqrt(nSpecies), Model=MODEL)
imp_cv_tss[imp_cv_tss$Var==ord,]
imp_cv_tss$upr<-imp_cv_tss$Imp + imp_cv_tss$se
imp_cv_tss$lwr<-imp_cv_tss$Imp - imp_cv_tss$se
assign(paste0('imp_cv_tss_', MODEL), imp_cv_tss)

#TSS pres
imp_pres_tss<-data.frame(Var=names(TSS_pres_Imp_mn), Imp=as.numeric(TSS_pres_Imp_mn), se=as.numeric(TSS_pres_Imp_sd)/sqrt(nSpecies), Model=MODEL)
imp_pres_tss[imp_pres_tss$Var==ord,]
imp_pres_tss$upr<-imp_pres_tss$Imp + imp_pres_tss$se
imp_pres_tss$lwr<-imp_pres_tss$Imp - imp_pres_tss$se
assign(paste0('imp_pres_tss_', MODEL), imp_pres_tss)

#TSS fut
imp_fut_tss<-data.frame(Var=names(TSS_fut_Imp_mn), Imp=as.numeric(TSS_fut_Imp_mn), se=as.numeric(TSS_fut_Imp_sd)/sqrt(nSpecies), Model=MODEL)
imp_fut_tss[imp_fut_tss$Var==ord,]
imp_fut_tss$upr<-imp_fut_tss$Imp + imp_fut_tss$se
imp_fut_tss$lwr<-imp_fut_tss$Imp - imp_fut_tss$se
assign(paste0('imp_fut_tss_', MODEL), imp_fut_tss)

load(paste0('Importances_',MODEL,'_AUC_R1.Rdata'))

#AUC cv
imp_cv_auc<-data.frame(Var=names(AUC_cv_Imp_mn), Imp=as.numeric(AUC_cv_Imp_mn), se=as.numeric(AUC_cv_Imp_sd)/sqrt(nSpecies), Model=MODEL)
imp_cv_auc[imp_cv_auc$Var==ord,]
imp_cv_auc$upr<-imp_cv_auc$Imp + imp_cv_auc$se
imp_cv_auc$lwr<-imp_cv_auc$Imp - imp_cv_auc$se
assign(paste0('imp_cv_auc_', MODEL), imp_cv_auc)

#AUC pres
imp_pres_auc<-data.frame(Var=names(AUC_pres_Imp_mn), Imp=as.numeric(AUC_pres_Imp_mn), se=as.numeric(AUC_pres_Imp_sd)/sqrt(nSpecies), Model=MODEL)
imp_pres_auc[imp_pres_auc$Var==ord,]
imp_pres_auc$upr<-imp_pres_auc$Imp + imp_pres_auc$se
imp_pres_auc$lwr<-imp_pres_auc$Imp - imp_pres_auc$se
assign(paste0('imp_pres_auc_', MODEL), imp_pres_auc)

#AUC fut
imp_fut_auc<-data.frame(Var=names(AUC_fut_Imp_mn), Imp=as.numeric(AUC_fut_Imp_mn), se=as.numeric(AUC_fut_Imp_sd)/sqrt(nSpecies), Model=MODEL)
imp_fut_auc[imp_fut_auc$Var==ord,]
imp_fut_auc$upr<-imp_fut_auc$Imp + imp_fut_auc$se
imp_fut_auc$lwr<-imp_fut_auc$Imp - imp_fut_auc$se
assign(paste0('imp_fut_auc_', MODEL), imp_fut_auc)

}

#MERGE MODELS
imp_cv_tss<-rbind(imp_cv_tss_RF, imp_cv_tss_GLM, imp_cv_tss_MaxEnt)
imp_cv_tss$Var <- factor(imp_cv_tss$Var, levels=ord); levels(imp_cv_tss$Var)<-labs
imp_cv_tss$Model<-factor(imp_cv_tss$Model, levels=c('GLM', 'MaxEnt', 'RF'))

imp_pres_tss<-rbind(imp_pres_tss_RF, imp_pres_tss_GLM, imp_pres_tss_MaxEnt)
imp_pres_tss$Var <- factor(imp_pres_tss$Var, levels=ord); levels(imp_pres_tss$Var)<-labs
imp_pres_tss$Model<-factor(imp_pres_tss$Model, levels=c('GLM', 'MaxEnt', 'RF'))

imp_fut_tss<-rbind(imp_fut_tss_RF, imp_fut_tss_GLM, imp_fut_tss_MaxEnt)
imp_fut_tss$Var <- factor(imp_fut_tss$Var, levels=ord); levels(imp_fut_tss$Var)<-labs
imp_fut_tss$Model<-factor(imp_fut_tss$Model, levels=c('GLM', 'MaxEnt', 'RF'))

imp_cv_auc<-rbind(imp_cv_auc_RF, imp_cv_auc_GLM, imp_cv_auc_MaxEnt)
imp_cv_auc$Var <- factor(imp_cv_auc$Var, levels=ord); levels(imp_cv_auc$Var)<-labs
imp_cv_auc$Model<-factor(imp_cv_auc$Model, levels=c('GLM', 'MaxEnt', 'RF'))

imp_pres_auc<-rbind(imp_pres_auc_RF, imp_pres_auc_GLM, imp_pres_auc_MaxEnt)
imp_pres_auc$Var <- factor(imp_pres_auc$Var, levels=ord); levels(imp_pres_auc$Var)<-labs
imp_pres_auc$Model<-factor(imp_pres_auc$Model, levels=c('GLM', 'MaxEnt', 'RF'))

imp_fut_auc<-rbind(imp_fut_auc_RF, imp_fut_auc_GLM, imp_fut_auc_MaxEnt)
imp_fut_auc$Var <- factor(imp_fut_auc$Var, levels=ord); levels(imp_fut_auc$Var)<-labs
imp_fut_auc$Model<-factor(imp_fut_auc$Model, levels=c('GLM', 'MaxEnt', 'RF'))


p1_tss<-ggplot(data=imp_cv_tss, aes(x=Var, y=Imp, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=lwr, ymax=upr), alpha=0.6, size=0.4, width=.4, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle('TSS cv') + xlab('') + ylab('% Variable importance') +
  scale_fill_brewer(palette="Paired") +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.y= element_text(colour='black'), legend.title=element_blank())

p2_tss<-ggplot(data=imp_pres_tss, aes(x=Var, y=Imp, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=lwr, ymax=upr), alpha=0.6, size=0.4, width=.4, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle('TSS present') + xlab('') + ylab('% Variable importance') +
  scale_fill_brewer(palette="Paired") +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.y= element_text(colour='black'), legend.title=element_blank())

p3_tss<-ggplot(data=imp_fut_tss, aes(x=Var, y=Imp, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=lwr, ymax=upr), alpha=0.6, size=0.4, width=.4, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle('TSS future') + xlab('') + ylab('% Variable importance') +
  scale_fill_brewer(palette="Paired") +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.y= element_text(colour='black'), legend.title=element_blank())
#
p1_auc<-ggplot(data=imp_cv_auc, aes(x=Var, y=Imp, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=lwr, ymax=upr), alpha=0.6, size=0.4, width=.4, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle('AUC cv') + xlab('') + ylab('% Variable importance') +
  scale_fill_brewer(palette="Paired") +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.y= element_text(colour='black'), legend.title=element_blank())

p2_auc<-ggplot(data=imp_pres_auc, aes(x=Var, y=Imp, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=lwr, ymax=upr), alpha=0.6, size=0.4, width=.4, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle('AUC present') + xlab('') + ylab('% Variable importance') +
  scale_fill_brewer(palette="Paired") +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.y= element_text(colour='black'), legend.title=element_blank())

p3_auc<-ggplot(data=imp_fut_auc, aes(x=Var, y=Imp, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=lwr, ymax=upr), alpha=0.6, size=0.4, width=.4, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle('AUC future') + xlab('') + ylab('% Variable importance') +
  scale_fill_brewer(palette="Paired") +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.y= element_text(colour='black'), legend.title=element_blank())


pdf('/Figs/Barplot_combined.pdf', width=13, height=8)
ggarrange(p1_tss, p2_tss, p3_tss, p1_auc, p2_auc, p3_auc,
	legend='top', common.legend=TRUE, labels=c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'), font.label=list(size=10, face='plain'),
	ncol=3, nrow=2)
	#ncol=3, nrow=2)
dev.off()



#EXP AND CONTRACTION AREA

library(ggplot)
library(ggpubr)

ord<-rev(c('SpeciesPrevalence_gExt','SamplePrevalence','nPresSampled', 'MESSland_ave','Bias', 'propNonEq','bufferP','TruePredictors','FalsePredictors'))
labs<-rev(c('Species Prevalence', 'Sample Prevalence', 'n Presences', 'Env. similarity', 'Env. gradient sampled', 'Niche Filling', 'Geographic Extent', 'Relevant Predictors', 'Irrelevant Predictors'))

models<-c('RF', 'MaxEnt', 'GLM')
for (i in 1:length(models)) {
MODEL=models[i]
nSpecies=50
load(paste0('Importances_',MODEL,'_R1.Rdata'))

#GAIN
imp_gain<-data.frame(Var=names(est_gain_per_tss_Imp_mn), Imp=as.numeric(est_gain_per_tss_Imp_mn), se=as.numeric(est_gain_per_tss_Imp_mn)/sqrt(nSpecies), Model=MODEL)
imp_gain[imp_gain$Var==ord,]
imp_gain$upr<-imp_gain$Imp + imp_gain$se
imp_gain$lwr<-imp_gain$Imp - imp_gain$se
assign(paste0('imp_gain_', MODEL), imp_gain)

#LOSS
imp_loss<-data.frame(Var=names(est_loss_per_tss_Imp_mn), Imp=as.numeric(est_loss_per_tss_Imp_mn), se=as.numeric(est_loss_per_tss_Imp_mn)/sqrt(nSpecies), Model=MODEL)
imp_loss[imp_loss$Var==ord,]
imp_loss$upr<-imp_loss$Imp + imp_loss$se
imp_loss$lwr<-imp_loss$Imp - imp_loss$se
assign(paste0('imp_loss_', MODEL), imp_loss)

}

#MERGE MODELS
imp_gain<-rbind(imp_gain_RF, imp_gain_GLM, imp_gain_MaxEnt)
imp_gain$Var <- factor(imp_gain$Var, levels=ord); levels(imp_gain$Var)<-labs
imp_gain$Model<-factor(imp_gain$Model, levels=c('GLM', 'MaxEnt', 'RF'))

imp_loss<-rbind(imp_loss_RF, imp_loss_GLM, imp_loss_MaxEnt)
imp_loss$Var <- factor(imp_loss$Var, levels=ord); levels(imp_loss$Var)<-labs
imp_loss$Model<-factor(imp_loss$Model, levels=c('GLM', 'MaxEnt', 'RF'))

p1<-ggplot(data=imp_gain, aes(x=Var, y=Imp, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=lwr, ymax=upr), alpha=0.6, size=0.4, width=.4, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle('Range expansion %') + xlab('') + ylab('% Variable importance') +
  scale_fill_brewer(palette="Paired") +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.y= element_text(colour='black'), legend.title=element_blank())

p2<-ggplot(data=imp_loss, aes(x=Var, y=Imp, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=lwr, ymax=upr), alpha=0.6, size=0.4, width=.4, position=position_dodge(0.9)) +
  coord_flip() +
  ggtitle('Range contraction %') + xlab('') + ylab('% Variable importance') +
  scale_fill_brewer(palette="Paired") +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.y= element_text(colour='black'), legend.title=element_blank())


pdf('/Figs/Fig. 6 Barplot Expansion Contraction.pdf', width=8, height=4)
ggarrange(p1, p2,
  legend='top', common.legend=TRUE, labels=c('(a)', '(b)'), font.label=list(size=10, face='plain'),
  ncol=2, nrow=1)
  #ncol=3, nrow=2)
dev.off()


