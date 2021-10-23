##########################################################################################################################################
#                                                                                                                                        #
# Author: XXXX XXXX                                                                                                                	     #
# e-mail: XXXXXX@XXXXX                                                                                                    				 #
# Reference:                                                                                                                             #
# XXXX et al. 2019. Assessing the reliability of species distribution projections in climate change research. XXXXXX                     #
#                                                                                                                                        #
##########################################################################################################################################

#Plot TSS vs AUC

RF<-read.csv('OUTPUT_Simulation_RF_R1.csv')
MX<-read.csv('OUTPUT_Simulation_MaxEnt_R1.csv')
GLM<-read.csv('OUTPUT_Simulation_GLM_R1.csv')

RF$n<-1:nrow(RF)
GLM$n<-1:nrow(GLM)
MX$n<-1:nrow(MX)

#clean
RF<-RF[complete.cases(RF[,c('SpID', 'TSS_cv', 'TSS_pres', 'TSS_fut', 'AUC_cv', 'AUC_pres', 'AUC_fut'),]),]
MX<-MX[complete.cases(MX[,c('SpID', 'TSS_cv', 'TSS_pres', 'TSS_fut', 'AUC_cv', 'AUC_pres', 'AUC_fut'),]),]
GLM<-GLM[complete.cases(GLM[,c('SpID', 'TSS_cv', 'TSS_pres', 'TSS_fut', 'AUC_cv', 'AUC_pres', 'AUC_fut'),]),]

set.seed(1)
rows<-sample(nrow(RF), 100000)
RF<-RF[rows,]
MX<-MX[rows,]
GLM<-GLM[rows,]

pc=19
cx=1
cxl=1.5
cxax=2
cxm=2

cl=adjustcolor('black', alpha=0.1)
tiff('Figs/AUC vs TSS.tif', height=1000, width=1100)
par(mfrow=c(3,3), mar=c(7, 7, 3, 2), las=1)
plot(GLM$TSS_cv, GLM$AUC_cv, cex.axis=cxax, cex.main=cxm, pch=pc, col=cl, bg=cl, cex=cx, main='GLM', xlab='', ylab='')
mtext('TSSpres', side=1, line=4, cex=cxl); mtext('AUCpres', side=2, line=5, cex=cxl, las=0); 
plot(GLM$TSS_pres, GLM$AUC_pres, cex.axis=cxax, cex.main=cxm, pch=pc, col=cl, bg=cl, cex=cx, main='GLM', xlab='', ylab='')
mtext('TSSfut', side=1, line=4, cex=cxl); mtext('AUCfut', side=2, line=5, cex=cxl, las=0); 
plot(GLM$TSS_fut, GLM$AUC_fut, cex.axis=cxax, cex.main=cxm, pch=pc, col=cl, bg=cl, cex=cx, main='GLM', xlab='', ylab='')
mtext('TSScv', side=1, line=4, cex=cxl); mtext('AUCcv', side=2, line=5, cex=cxl, las=0); 
plot(MX$TSS_cv, MX$AUC_cv, cex.axis=cxax, cex.main=cxm, pch=pc, col=cl, bg=cl, cex=cx, main='MaxEnt', xlab='', ylab='')
mtext('TSScv', side=1, line=4, cex=cxl); mtext('AUCcv', side=2, line=5, cex=cxl, las=0); 
plot(MX$TSS_pres, MX$AUC_pres, cex.axis=cxax, cex.main=cxm, pch=pc, col=cl, bg=cl, cex=cx, main='MaxEnt', xlab='', ylab='')
mtext('TSSpres', side=1, line=4, cex=cxl); mtext('AUCpres', side=2, line=5, cex=cxl, las=0); 
plot(MX$TSS_fut, MX$AUC_fut, cex.axis=cxax, cex.main=cxm, pch=pc, col=cl, bg=cl, cex=cx, main='MaxEnt', xlab='', ylab='')
mtext('TSSfut', side=1, line=4, cex=cxl); mtext('AUCfut', side=2, line=5, cex=cxl, las=0); 
plot(RF$TSS_cv, RF$AUC_cv, cex.axis=cxax, cex.main=cxm, pch=pc, col=cl, bg=cl, cex=cx, main='Random Forest', xlab='', ylab='')
mtext('TSScv', side=1, line=4, cex=cxl); mtext('AUCcv', side=2, line=5, cex=cxl, las=0); 
plot(RF$TSS_pres, RF$AUC_pres, cex.axis=cxax, cex.main=cxm, pch=pc, col=cl, bg=cl, cex=cx, main='Random Forest', xlab='', ylab='')
mtext('TSSpres', side=1, line=4, cex=cxl); mtext('AUCpres', side=2, line=5, cex=cxl, las=0); 
plot(RF$TSS_fut, RF$AUC_fut, cex.axis=cxax, cex.main=cxm, pch=pc, col=cl, bg=cl, cex=cx, main='Random Forest', xlab='', ylab='')
mtext('TSSfut', side=1, line=4, cex=cxl); mtext('AUCfut', side=2, line=5, cex=cxl, las=0); 
dev.off()

