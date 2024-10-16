#this script is for running enrichr with DEP data and visualize them 
#1. Load DEP file, subset DEGs by a defined p-value threshold, run EnrichR
#2. Query STRING for significant TransFac TFs + associated DEGs
#3. Run MCode on STRING network, run Enrichr on top clusters

library(openxlsx)
library(enrichR)
library(glue)
library(rlist)
library(STRINGdb)
library(dplyr)
library(forcats)

setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/GO_ontology/")
setwd("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/DEP&VOLCANO/Qns_only")
#loading DEP
DEPall = read.xlsx("DEP.onlyQueens(15subj).withoutSANDY.xlsx", sheet = 2)
upDEP = read.xlsx("DEP.onlyQueens(15subj).withoutSANDY.xlsx", sheet = 4)
downDEP = read.xlsx("DEP.onlyQueens(15subj).withoutSANDY.xlsx", sheet = 3)
#loading list of unique proteins from MS and Qns cohort
MS_unique_controls = read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/data/Unique.Protein.MS.xlsx", sheet = 1)
colnames(MS_unique_controls) = MS_unique_controls[3,]
MS_unique_controls = MS_unique_controls[4:nrow(MS_unique_controls),]
#loading unique proteins that are only present in Cannabis but not in controls
MS_unique_cannabis = read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/data/Unique.Protein.MS.xlsx", sheet = 2)
colnames(MS_unique_cannabis) = MS_unique_cannabis[3,]
MS_unique_cannabis = MS_unique_cannabis[4:nrow(MS_unique_cannabis),]
#loading unique proteins that are only present in Controls but not in cannabis in Queens
Qns_unique_control = read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/data/Unique.Protein.Queens.xlsx", sheet = 1)
colnames(Qns_unique_control) = Qns_unique_control[3,]
Qns_unique_control = Qns_unique_control[4:nrow(Qns_unique_control),]
#loading unique proteins that are only present in Cannabis but not in controls in Queens
Qns_unique_cannabis = read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/data/Unique.Protein.Queens.xlsx", sheet = 2)
colnames(Qns_unique_cannabis) = Qns_unique_cannabis[1,]
Qns_unique_cannabis = Qns_unique_cannabis[2:nrow(Qns_unique_cannabis),]

#loading DEP with controlling covariates
DEPall <- read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/DEP&VOLCANO/bothsex.DEP.nostress50human.xlsx", sheet = 3)
upDEP <- read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/DEP&VOLCANO/bothsex.DEP.nostress50human.xlsx", sheet = 4)
downDEP <- read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/DEP&VOLCANO/bothsex.DEP.nostress50human.xlsx", sheet = 5)
#for male
maleDEP <- read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/Male/proteomics/data/maleDEP.human.xlsx", sheet= 3)
maleup <- read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/Male/proteomics/data/maleDEP.human.xlsx", sheet= 4)
maledown <- read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/Male/proteomics/data/maleDEP.human.xlsx", sheet= 5)
#for female
femaleDEP <- read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/Female/proteomics/DEP/femaleDEP.human.xlsx", sheet = 3)
female_up <- read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/Female/proteomics/DEP/femaleDEP.human.xlsx", sheet = 4)
female_down <- read.xlsx("/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/Female/proteomics/DEP/femaleDEP.human.xlsx", sheet = 5)

#customize variables
pval_threshold = 0.05
symbol_column = 'Accession'
species = 'human'
####never get too much worried on setting working directory or the version, that will be changed according to the files need to process :)
#setwd("/Users/teestanaskar/Dropbox/Anissa=Randy=Teesta/Placenta_Project/RNA-seq/placenta/Results/deg_analyses/human/Teesta.Jan2023.WithoutSANDY.93subj/WithCovariatescontrolled/")
folder = "/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/GO_ontology/Queens_DEP/"
#for male 
folder = "/Users/teestanaskar/Dropbox/Teesta/Placenta/Human.Placenta/bothsex/Proteomics/GO_ontology/"
#for female
folder = "/Users/teestanaskar/Dropbox/Teesta/Placenta/Rat_placenta/Proteomics/GO_ontology/dep/enrichr_outputs/"
#making a dataframe only for genes those are above pvalue threshhold
#Enrichr
dbs = c('GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023',
        'KEGG_2021_Human', 'DisGeNET', 'TRANSFAC_and_JASPAR_PWMs')
output = enrichr(MS_unique_cannabis$Accession, databases=dbs)
output = enrichr(Qns_unique_cannabis$Accession, databases=dbs)
output = enrichr(DEP[DEP$pvalue < 0.05,]$Accession, databases=dbs)
output = enrichr(MS_DEP[MS_DEP$P.Value< 0.05,]$Accession, databases = dbs)
output = enrichr(DEPall$Accession, databases=dbs)
#celltissuetypes = output$`Descartes Cell Types and Tissue 2021`
#Save Enrichr output as Excel file
#for unique 
output = list.prepend(output, Qns_unique_cannabis$Accession)
filename = "human_Qns_unique_in_cannabis.ontologies.xlsx"
write.xlsx(output, glue(folder, filename))

output = list.prepend(output, DEP[DEP$pvalue<0.05,]$Gene.Names)
filename = "ratDEP.bothsex.all.Accession.ontologies.xlsx"
write.xlsx(output, glue(folder, filename))
##for DEP queens
output = list.prepend(output, DEPall$Accession)
filename = "human_DEP_queens.ontologies.xlsx"
write.xlsx(output, glue(folder, filename))

#visualization of enrichmentn result
#worked the ggplot bar plot
setwd("/Users/teestanaskar/Dropbox/Anissa=Randy=Teesta/Placenta_Project/RNA-seq/placenta/Results/deg_analyses/human/Teesta.Jan2023.WithoutSANDY.93subj/WithCovariatescontrolled/Ontology_for_SigDEGs.total.UP.DOWN")
BP.down <- read.xlsx("ontology/Ontology.DOWN.DEG.xlsx", sheet = 2)
BP.up <- read.xlsx("UP.sig.COVariates.cntrld.DEG.93placenta.xlsx", sheet = 2)
Immune <- BP.down[BP.down$Category == "Immune" & BP.down$P.value< 0.05,]

Immune <- BP.up[BP.up$Category== "Immune" & BP.up$P.value< 0.05,]
BP.down$`-log10P` = -log10(BP.down$P.value)
BP.down[1:10,] %>% mutate(Term = fct_reorder(Term, `-log10P`)) %>% ggplot(aes(x=Term, `-log10P`)) + geom_bar(stat="identity", fill="red", alpha=1, width=.7) +
  coord_flip() +
  xlab("") +
  theme_bw()


BP <- output$GO_Biological_Process_2018
#sorting
DisGeNET <- DisGeNET[order(DisGeNET$Adjusted.P.value), ]
DisGeNET <- DisGeNET[order(DisGeNET$P.value),]
#visualizing
DisGeNET$minuslog.Padj <- -log10(DisGeNET$Adjusted.P.value)
DisGeNET$minuslog.Padj <- -log10(DisGeNET$Adjusted.P.value)
DisGeNET[1:20,] %>% mutate(Term = fct_reorder(Term, minuslog.Padj)) %>% ggplot(aes(x=Term, minuslog.Padj)) + geom_bar(stat="identity", fill="black", alpha=1, width=.7) +
  coord_flip() +
  xlab("") +
  theme_bw()

DEP_bothsex_ratPlacenta_DisGeNET$minuslog.Padj <- -log10(DEP_bothsex_ratPlacenta_DisGeNET$Adjusted.P.value)
DEP_bothsex_ratPlacenta_DisGeNET <- DEP_bothsex_ratPlacenta_DisGeNET[1:20,] %>% mutate(Term = fct_reorder(Term, minuslog.Padj)) %>% ggplot(aes(x=Term, minuslog.Padj)) + geom_bar(stat="identity", fill="black", alpha=1, width=.7) +
  coord_flip() +
  xlab("") +
  theme_bw()

Male.TFENCODE[1:20,] %>% mutate(Term = fct_reorder(Term, Odds.Ratio)) %>% ggplot(aes(x=Term, Odds.Ratio)) + geom_bar(stat="identity", fill="black", alpha=1, width=.7) +
  coord_flip() +
  xlab("") +
  theme_bw()
tiff("GO_BP.image.tiff, res= 1200, width =3000, height = 5500")
print(GO_BP.image)
dev.off()
RStudioGD 
2 

# string_query = c()
# for (i in 1:nrow(transfac_sig)) {
#   tf = strsplit(output$TRANSFAC_and_JASPAR_PWMs[i,1]," ", fixed=-T)[[1]][1]
#   tf_genes = strsplit(transfac_sig[i,ncol(transfac_sig)],";")[[1]]
#   string_query = c(string_query, tf, tf_genes)
# }

#STRINGdb
#species: Human: 9606, Rat: 10116, Mouse: 10090
if(species == 'rat'){
  string_species = 10116 
} else if(species == 'human'){
  string_species = 9606
} else if(species == 'mouse'){
  string_species = 10090
}

string_query = genes[,1]
score_thresh = 700
folder = Users/teestanaskar/Desktop/HurdLab.Data/prenatal_cannabis_use/Rat_placenta_proteomics_Enrichment/bothsex
string_db <- STRINGdb$new( version="11", species=string_species, score_threshold=score_thresh, input_directory="")
string_map <- string_db$map(as.data.frame(string_query), 'string_query', removeUnmappedRows = T)
graph <- string_db$get_interactions(string_map$STRING_id)
df = left_join(graph, string_map, by=c("from" = "STRING_id"))
df = left_join(df, string_map, by=c("to" = "STRING_id"))
write.csv(df, glue(/Users/teestanaskar/Desktop/HurdLab.Data/prenatal_cannabis_use/Rat_placenta_proteomics_Enrichment/bothsex.DEG, 'string/', 'string_network_','allDEP.bothsex_','score',score_thresh,'.csv'))

