#this script is for human placental differential expression of proteins analysis by controlling covariates
setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/data/")
library(limma)
df.prot <- read.xlsx("normailzed_onlyQueens(withoutSANDY).xlsx", sheet=4) #for queens
df.prot <- read.xlsx("normalized_onlyMS_withoutSANDY_bothsexcombined.xlsx", sheet = 4)
data = df.prot[,4:ncol(df.prot)]
rownames(data) = df.prot$Accession
#load the meta data
meta <- read.xlsx("/Users/teestanaskar/Dropbox/Anissa=Randy=Teesta/Placenta_Project/HUMAN.PLACENTA.METADATA/WorkingMetadata.updatedbyTeesta.consolidatedinfo.Anissa.Greg.Placenta.inventory.xlsx", sheet = 4)
###removed samples whos group status is not matches among Anissa and Greg
#meta = meta[! meta$Placenta_Seq_ID %in% c('S_282M', 'S_29', 'S_67', 'S_137'),]
sub.col <- meta[meta$Site=="Qns",]$YaleID.Proteomics # when queens only
sub.col <- meta[meta$Site=="MS",]$YaleID.Proteomics #when MS only
dat <- data[, sub.col]
rownames(dat) <- rownames(data)
dat.log = log2(dat+1)
#there were error saying about characters peresnt on the data , so converting it to numeric first
dat[] <- lapply(dat, as.numeric)
dat.log = log2(dat+1)
#dat.log = na.omit(dat.log)
#boxplot(dat,las=2,main= "dat")
#dat = equalMedianNormalization(dat)
boxplot(dat.log,las=2,main= "dat")
meta=meta[meta$Site=="Qns",][match(colnames(dat.log),meta[meta$Site=="Qns",]$YaleID.Proteomics),]
meta=meta[meta$Site=="MS",][match(colnames(dat.log),meta[meta$Site=="MS",]$YaleID.Proteomics),]
#meta=meta[match(colnames(normalized_data),meta$YaleID.Proteomics),]
rownames(meta) <- meta$YaleID.Proteomics
all(rownames(meta) %in% colnames(dat.log)) #TRUE
#checking whether the order are same
all(rownames(meta) == colnames(dat.log)) #TRUE
meta = meta[,c("Group", "CSEX", "smoking", "DrinkD", "vag_infec", "GA", "MOB_AGE")]
group = as.factor(meta$Group)
#site = as.factor(meta$Site) keep it off while doing for only one site
CSEX = as.factor(meta$CSEX)
smoke = as.factor(meta$smoking)
drink = as.factor(meta$DrinkD)
vaginfect = as.factor(meta$vag_infec)
GA= as.numeric(meta$GA)
age = as.numeric(meta$MOB_AGE)
#creating a design matrix for measuring cannabis over control on controlling mother's age, child sex, gestational age, smoke, drink and vaginal infection              
design = model.matrix(~0 +group+age+CSEX+GA+smoke+drink+vaginfect)
colnames(design) = gsub("group","",colnames(design))
#renaming
colnames(design) <- c("Cannabis", "Control", "MOB_AGE", "CSEX", "GA", "smoke", "drink", "vaginfect")
#making contrast
contr.matrix <- makeContrasts(ControlvsCannabis = Control-Cannabis, 
                              levels = colnames(design))
contr.matrix <- makeContrasts(Cannabis - Control, levels = design)
fit1 <- lmFit(dat.log, design)
fit2 <- contrasts.fit(fit1, contr.matrix)
ebFit <- eBayes(fit2)
limma.results <- topTable(ebFit, coef=1, n= Inf)
limma.results$Accession = rownames(limma.results)
limma.results$abundance = dat[limma.results$Accession,]$Freq
dt <- decideTests(ebFit)
summary(dt)
## Add direction of log fold change relative to control
limma.results$direction <- ifelse(limma.results$logFC > 0,
                                  "up", "down") %>%
  as.factor()

## Add significance thresholds
limma.results$significance <- ifelse(limma.results$P.Value < 0.05,
                                     "sig", "not.sig") %>%
  as.factor()
limma.results_withdata <- inner_join(limma.results, df.prot, by = "Accession")
write.xlsx(limma.results, "../DEP&VOLCANO/DEP.MSonly(35subj).covariatescontrolled.xlsx")
write.csv(limma.results_withdata, "../DEP&VOLCANO/bothsexDEP.nostress50humanplacenta.csv")
plotMD(ebFit, column=1, status=dt[,1]) 
volcanoplot(limma.results)
limma::plotMA(limma.results)
