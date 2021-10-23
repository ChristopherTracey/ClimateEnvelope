##########################################################################################################################################
#                                                                                                                                        #
# Author: XXXX XXXX                                                                                                                	     #
# e-mail: XXXXXX@XXXXX                                                                                                    				 #
# Reference:                                                                                                                             #
# XXXX et al. 2019. Assessing the reliability of species distribution projections in climate change research. XXXXXX                     #
#                                                                                                                                        #
##########################################################################################################################################

#This code combines all array job outputs in one single file

library(plyr)

MODEL='GLM'
files<-dir('/OutArrayJob/', recursive=FALSE, full.names=TRUE)
num<-grep(MODEL, files)
files<-files[num]

lista<-list(0)
for (i in 1:length(files)){
lista[[i]]<-read.csv(as.character(files[i]), header=TRUE)
}

out<-rbind.fill(lista)

write.csv(out, paste0('/vol/milkunB/lsantini/Vspecies/OUTPUT_Simulation_',MODEL,'_R1.csv'), row.names=FALSE)

