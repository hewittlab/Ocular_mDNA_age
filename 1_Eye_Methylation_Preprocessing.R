#############################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("minfi")
#biocLite("IlluminaHumanMethylation450kmanifest")
#############################################################

library(minfi)
rgset<- read.450k.exp("~/Desktop/Research/AMD/Methylation/EyeBank/Data/raw/Eye Methylation Project/")
rgset
pData(rgset)
head(sampleNames(rgset))
beta.table<-getBeta(rgset)
colnames(beta.table)
colnames(beta.table) <- c("3631B","3675O","3631C","3675R","3631O","3677B","3631R","3677C","3675B","3677O","3675C","3677R","3684B","3685O","3684C","3685R","3684O","3689B","3684R","3689C","3685B","3689O","3685C","3689R","3682C","3682O","3682R","3701R","3701B")
write.csv(beta.table,"eyebetas.csv", row.names = TRUE)
