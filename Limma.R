setwd("/home/dm/Desktop/Recent_work/Bioinformatics/Scleroderma/SSc/GSE104174")

# 2. Load libraries for script one
library(RCurl)
library(GEOquery)
library(limma)

url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20966/matrix/"



gse <- getGEO(filename = "GSE104174_series_matrix.txt.gz",destdir=getwd())
d <- factor(c(rep('CTRL', 7),rep('T2D',10)))


#d <- factor(c('CTRL','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','CTRL','T2D','CTRL','T2D','T2D','CTRL','T2D','CTRL','T2D','CTRL','T2D','CTRL','CTRL','T2D','T2D'))


mod <- model.matrix(~0+d)
# print(mod)
# dCD dCTRL dCVCTRL dCVID
# 1    0     0       0     1
# 2    0     0       1     0

fit_1 <- lmFit(gse, mod)
# #print(fit_1)
contr <- makeContrasts(dCTRL-dT2D,levels = mod)
#contr <- makeContrasts(dCD-dCTRL, levels = mod)
fit_2 <- contrasts.fit(fit = fit_1, contrasts = contr)
fit_3 <- eBayes(fit_2)


#[1] 156  32
table_result <- topTable(fit_3, coef=1,n=Inf,adjust="BH", sort.by = "logFC")
