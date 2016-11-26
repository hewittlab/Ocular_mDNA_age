library(WGCNA)
library(sqldf)
library(impute)
library(lattice)

setwd("~/Desktop/Research/AMD/Methylation/EyeBank/Data/analysis/")

dataout <- read.csv("EyeAge_Output.csv")
colnames(dataout)
dataout <- dataout[order(dataout$SampleID),] 
signif(dataout$DNAmAge,2)

datSample<-read.csv("sampleidannotat.csv")
colnames(datSample)
datSample <- datSample[order(datSample$id),] 
DNAmAge=dataout$DNAmAge
medianAbsDev=function(x,y) median(abs(x-y),na.rm=TRUE)
medianAbsDev1=signif(medianAbsDev(DNAmAge, datSample$Age),2)
par(mfrow=c(1,1))

#### Figure 1a
#verboseScatterplot(datSample$Age,DNAmAge,xlab="Chronological Age",ylab="DNAm Age",ylim=c(10,100), pch =c(1,2,3,4)[as.numeric(datSample$SampleLabel)], cex = 1.3, cex.axis = 1, cex.lab = 1, cex.main = 1, col= rainbow(10)[round(datSample$Age/max(datSample$Age)*10)], main=paste("All, err=", medianAbsDev1) ) ;abline(0,1) 
verboseScatterplot(datSample$Age,DNAmAge,xlab="Chronological Age",ylab="DNAm Age",xlim=c(10,80),ylim=c(10,80), pch =c(1,2,3,4)[as.numeric(datSample$SampleLabel)], cex = 1.8, cex.axis = 1, cex.lab = 1,cex.main=1,col.main="white",col= c("black"));abline(0,1) 
#verboseScatterplot(datSample$Age,DNAmAge,xlab="Chronological Age",ylab="DNAm Age",ylim=c(10,100), pch =c(15,16,17,18)[as.numeric(datSample$SampleLabel)], cex = 1.3, cex.axis = 1, cex.lab = 1, cex.main = 1, col= c("black"), main=paste("All, err=", medianAbsDev1) ) ;abline(0,1) 

legend("topleft", legend= c("Blood","Neurosensory Retina","RPE/Choriod","Optic Nerve"),pch=c(1,2,3,4),bty = "n" ,cex = 1)

#### Figure 1b

AgeaccelerationDiff = DNAmAge - datSample$Age
samlabel <- factor(c(datSample$SampleLabel) )
samlabel <- ordered(datSample$SampleLabel, levels=c("Blood","Optic Nerve","RPE/Choroid","Neurosensory Retina") )

verboseBarplot(AgeaccelerationDiff,samlabel, xlab= "Tissues", ylab = "Difference in epigenetic and chronological age(years)",ylim=c(-60,60), cex = 1, cex.axis = 1, cex.lab = 1, cex.main =1,col.main="white", numberStandardErrors = 1)

#### Figure 1c

#verboseBarplot(AgeaccelerationDiff, datSample$CauseofDeath,main = "Cause of Death", xlab ="Cause of Death", ylab = "Age Acceleration Diff",ylim=c(-70,70), cex =0.7, cex.axis =1,cex.lab =1, cex.main =1, numberStandardErrors = 1)
verboseScatterplot(datSample$PreservationTimeInterval,AgeaccelerationDiff,xlab="Post Mortem Interval",xlim=c(2,12),ylim=c(-60,20), pch =c(1,2,3,4),cex = 1.8, cex.axis =1,cex.lab =1,cex.main =1,main ="C",col.main="white")
