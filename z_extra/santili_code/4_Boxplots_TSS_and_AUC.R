##########################################################################################################################################
#                                                                                                                                        #
# Author: XXXX XXXX                                                                                                                      #
# e-mail: XXXXXX@XXXXX                                                                                                             #
# Reference:                                                                                                                             #
# XXXX et al. 2019. Assessing the reliability of species distribution projections in climate change research. XXXXXX                     #
#                                                                                                                                        #
##########################################################################################################################################

#Plot distribution of TSS and AUC values

library(ggplot2)
library(ggpubr)

models<-c('RF', 'MaxEnt', 'GLM')
for (i in 1:length(models)) {
MODEL=models[i]

data<-read.csv(paste0('OUTPUT_Simulation_',MODEL,'_R1.csv'))

D1=data[data$nPresSampled>200 & data$propNonEq==100 & data$Bias==100 & data$FalsePredictors==0,]
D2=data[data$nPresSampled<50 & data$propNonEq<50 & data$Bias<50 & data$TruePredictors==0,]

nrow(D1)
nrow(D2)

#TSS
D1tss<-D1
D2tss<-D2

TSS_cv<-data.frame(Var='TSS_cv', Value=D1tss$TSS_cv)
TSS_pres<-data.frame(Var='TSS_pres', Value=D1tss$TSS_pres)
TSS_fut<-data.frame(Var='TSS_fut', Value=D1tss$TSS_fut)
TSS_HabitatLoss<-data.frame(Var='TSS_HabitatLoss', Value=D1tss$TSS_HabitatLoss)
TSS_HabitatGain<-data.frame(Var='TSS_HabitatGain', Value=D1tss$TSS_HabitatGain)

D1tss<-rbind(TSS_cv, TSS_pres, TSS_fut, TSS_HabitatLoss, TSS_HabitatGain)

TSS_cv<-data.frame(Var='TSS_cv', Value=D2tss$TSS_cv)
TSS_pres<-data.frame(Var='TSS_pres', Value=D2tss$TSS_pres)
TSS_fut<-data.frame(Var='TSS_fut', Value=D2tss$TSS_fut)
TSS_HabitatLoss<-data.frame(Var='TSS_HabitatLoss', Value=D2tss$TSS_HabitatLoss)
TSS_HabitatGain<-data.frame(Var='TSS_HabitatGain', Value=D2tss$TSS_HabitatGain)

D2tss<-rbind(TSS_cv, TSS_pres, TSS_fut, TSS_HabitatLoss, TSS_HabitatGain)

#AUC
D1auc<-D1
D2auc<-D2

AUC_cv<-data.frame(Var='AUC_cv', Value=D1auc$AUC_cv)
AUC_pres<-data.frame(Var='AUC_pres', Value=D1auc$AUC_pres)
AUC_fut<-data.frame(Var='AUC_fut', Value=D1auc$AUC_fut)
AUC_HabitatLoss<-data.frame(Var='AUC_HabitatLoss', Value=D1auc$AUC_HabitatLoss)
AUC_HabitatGain<-data.frame(Var='AUC_HabitatGain', Value=D1auc$AUC_HabitatGain)

D1auc<-rbind(AUC_cv, AUC_pres, AUC_fut, AUC_HabitatLoss, AUC_HabitatGain)

AUC_cv<-data.frame(Var='AUC_cv', Value=D2auc$AUC_cv)
AUC_pres<-data.frame(Var='AUC_pres', Value=D2auc$AUC_pres)
AUC_fut<-data.frame(Var='AUC_fut', Value=D2auc$AUC_fut)
AUC_HabitatLoss<-data.frame(Var='AUC_HabitatLoss', Value=D2auc$AUC_HabitatLoss)
AUC_HabitatGain<-data.frame(Var='AUC_HabitatGain', Value=D2auc$AUC_HabitatGain)

D2auc<-rbind(AUC_cv, AUC_pres, AUC_fut, AUC_HabitatLoss, AUC_HabitatGain)

D1tss$Model<-MODEL
assign(paste0('D1tss_',MODEL), D1tss)

D2tss$Model<-MODEL
assign(paste0('D2tss_',MODEL), D2tss)

D1auc$Model<-MODEL
assign(paste0('D1auc_',MODEL), D1auc)

D2auc$Model<-MODEL
assign(paste0('D2auc_',MODEL), D2auc)

}

#COMBINE
D1tss<-rbind(D1tss_RF, D1tss_GLM, D1tss_MaxEnt)
D2tss<-rbind(D2tss_RF, D2tss_GLM, D2tss_MaxEnt)

D1auc<-rbind(D1auc_RF, D1auc_GLM, D1auc_MaxEnt)
D2auc<-rbind(D2auc_RF, D2auc_GLM, D2auc_MaxEnt)

#PLOT
levels(D1tss$Var) <- c('Estimated', 'Present', 'Future', 'Contraction', 'Expansion')
levels(D2tss$Var) <- c('Estimated', 'Present', 'Future', 'Contraction', 'Expansion')
levels(D1auc$Var) <- c('Estimated', 'Present', 'Future', 'Contraction', 'Expansion')
levels(D2auc$Var) <- c('Estimated', 'Present', 'Future', 'Contraction', 'Expansion')

D1tss$Model<-as.factor(D1tss$Model)
D2tss$Model<-as.factor(D2tss$Model)
D1auc$Model<-as.factor(D1auc$Model)
D2auc$Model<-as.factor(D2auc$Model)

p1tss<-ggplot(D1tss, aes(x=Var, y=Value, fill=Model)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle('Optimal conditions') + xlab('') + ylab('TSS') +
  coord_cartesian(ylim = c(-0.2, 1)) +
  #ylim(c(-0.2,1)) +
  #scale_fill_manual(values=cols) +
  scale_fill_brewer(palette='Paired') +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), legend.title=element_blank()) +
  geom_hline(yintercept=0, linetype=2) +
  geom_hline(yintercept=0.5, linetype=3)

p2tss<-ggplot(D2tss, aes(x=Var, y=Value, fill=Model)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle('Non-optimal conditions') + xlab('') + ylab('TSS') +
  coord_cartesian(ylim = c(-0.2, 1)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), legend.title=element_blank()) +
  geom_hline(yintercept=0, linetype=2) +
  geom_hline(yintercept=0.5, linetype=3)


p1auc<-ggplot(D1auc, aes(x=Var, y=Value, fill=Model)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle('Optimal conditions') + xlab('') + ylab('AUC') +
  coord_cartesian(ylim = c(0.45, 1)) +
  scale_fill_brewer(palette='Paired') +
  scale_y_continuous(breaks=seq(0.5, 1, 0.1)) +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), legend.title=element_blank()) +
  geom_hline(yintercept=0.5, linetype=2) +
  geom_hline(yintercept=0.7, linetype=3)

p2auc<-ggplot(D2auc, aes(x=Var, y=Value, fill=Model)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle('Non-optimal conditions') + xlab('') + ylab('AUC') +
  coord_cartesian(ylim = c(0.45, 1)) +
  scale_fill_brewer(palette='Paired') +
  scale_y_continuous(breaks=seq(0.5, 1, 0.1)) +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), legend.title=element_blank()) +
  geom_hline(yintercept=0.5, linetype=2) +
  geom_hline(yintercept=0.7, linetype=3)


pdf('~/Documents/Progetti/VirtualSDMclimateProjection/R1_Outputs/Figs/Boxplots_combined.pdf', width=10, height=8)
ggarrange(p1auc, p2auc, p1tss, p2tss, 
  legend='top', common.legend=TRUE, labels=c('(a)', '(b)', '(c)', '(d)'), font.label=list(size=10, face='plain'),
  ncol=2, nrow=2)
dev.off()


#############################################
############FIGS WITH NO SCENARIO ###########
#############################################

library(ggplot2)
library(ggpubr)

models<-c('RF', 'MaxEnt', 'GLM')
for (i in 1:length(models)) {
MODEL=models[i]

data<-read.csv(paste0('OUTPUT_Simulation_',MODEL,'_R1.csv'))

#TSS
Dtss<-data

TSS_cv<-data.frame(Var='TSS_cv', Value=Dtss$TSS_cv)
TSS_pres<-data.frame(Var='TSS_pres', Value=Dtss$TSS_pres)
TSS_fut<-data.frame(Var='TSS_fut', Value=Dtss$TSS_fut)
TSS_HabitatLoss<-data.frame(Var='TSS_HabitatLoss', Value=Dtss$TSS_HabitatLoss)
TSS_HabitatGain<-data.frame(Var='TSS_HabitatGain', Value=Dtss$TSS_HabitatGain)

Dtss<-rbind(TSS_cv, TSS_pres, TSS_fut, TSS_HabitatLoss, TSS_HabitatGain)


#AUC
Dauc<-data

AUC_cv<-data.frame(Var='AUC_cv', Value=Dauc$AUC_cv)
AUC_pres<-data.frame(Var='AUC_pres', Value=Dauc$AUC_pres)
AUC_fut<-data.frame(Var='AUC_fut', Value=Dauc$AUC_fut)
AUC_HabitatLoss<-data.frame(Var='AUC_HabitatLoss', Value=Dauc$AUC_HabitatLoss)
AUC_HabitatGain<-data.frame(Var='AUC_HabitatGain', Value=Dauc$AUC_HabitatGain)

Dauc<-rbind(AUC_cv, AUC_pres, AUC_fut, AUC_HabitatLoss, AUC_HabitatGain)

Dtss$Model<-MODEL
assign(paste0('Dtss_',MODEL), Dtss)

Dauc$Model<-MODEL
assign(paste0('Dauc_',MODEL), Dauc)

}

#COMBINE
Dtss<-rbind(Dtss_RF, Dtss_GLM, Dtss_MaxEnt)
Dauc<-rbind(Dauc_RF, Dauc_GLM, Dauc_MaxEnt)

#PLOT
levels(Dtss$Var) <- c('Estimated', 'Present', 'Future', 'Contraction', 'Expansion')
levels(Dauc$Var) <- c('Estimated', 'Present', 'Future', 'Contraction', 'Expansion')

Dtss$Model<-as.factor(Dtss$Model)
Dauc$Model<-as.factor(Dauc$Model)



ptss<-ggplot(Dtss, aes(x=Var, y=Value, fill=Model)) +
  geom_boxplot(outlier.shape=NA) +
  xlab('') + ylab('TSS') +
  coord_cartesian(ylim = c(-0.2, 1)) +
  scale_fill_brewer(palette='Paired') +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), legend.title=element_blank()) +
  geom_hline(yintercept=0, linetype=2) +
  geom_hline(yintercept=0.5, linetype=3)

pauc<-ggplot(Dauc, aes(x=Var, y=Value, fill=Model)) +
  geom_boxplot(outlier.shape=NA) +
  xlab('') + ylab('AUC') +
  coord_cartesian(ylim = c(0.45, 1)) +
  scale_fill_brewer(palette='Paired') +
  scale_y_continuous(breaks=seq(0.5, 1, 0.1)) +
  theme_classic() +
  theme(plot.title=element_text(size=10, hjust=0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=10), legend.title=element_blank()) +
  geom_hline(yintercept=0.5, linetype=2) +
  geom_hline(yintercept=0.7, linetype=3)


pdf('~/Documents/Progetti/VirtualSDMclimateProjection/R1_Outputs/Figs/Boxplots_combined_noScenarios.pdf', width=10, height=4)
ggarrange(pauc, ptss,
  legend='top', common.legend=TRUE, labels=c('(a)', '(b)'), font.label=list(size=10, face='plain'),
  ncol=2, nrow=1)
dev.off()
