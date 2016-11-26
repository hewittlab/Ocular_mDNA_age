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

# Sample Mapping
# Sentrix_ID_Position		ID
# 6042308131_R01C01	=	3631B
# 6042308131_R01C02	=	3675O
# 6042308131_R02C01	=	3631C
# 6042308131_R02C02	=	3675R
# 6042308131_R03C01	=	3631O
# 6042308131_R03C02	=	3677B
# 6042308131_R04C01	=	3631R
# 6042308131_R04C02	=	3677C
# 6042308131_R05C01	=	3675B
# 6042308131_R05C02	=	3677O
# 6042308131_R06C01	=	3675C
# 6042308131_R06C02	=	3677R
# 6042308133_R01C01	=	3684B
# 6042308133_R01C02	=	3685O
# 6042308133_R02C01	=	3684C
# 6042308133_R02C02	=	3685R
# 6042308133_R03C01	=	3684O
# 6042308133_R03C02	=	3689B
# 6042308133_R04C01	=	3684R
# 6042308133_R04C02	=	3689C
# 6042308133_R05C01	=	3685B
# 6042308133_R05C02	=	3689O
# 6042308133_R06C01	=	3685C
# 6042308133_R06C02	=	3689R
# 6042308134_R01C01	=	3682C
# 6042308134_R02C01	=	3682O
# 6042308134_R03C01	=	3682R
# 6042308134_R04C01	=	3701R
# 6042308134_R05C01	=	3701B

colnames(beta.table) <- c("3631B","3675O","3631C","3675R","3631O","3677B","3631R","3677C","3675B","3677O","3675C","3677R","3684B","3685O","3684C","3685R","3684O","3689B","3684R","3689C","3685B","3689O","3685C","3689R","3682C","3682O","3682R","3701R","3701B")
write.csv(beta.table,"eyebetas.csv", row.names = TRUE)
