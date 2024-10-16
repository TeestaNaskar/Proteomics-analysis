####this script is for visualization
library(EnhancedVolcano)
#customizing colors
keyvals <- ifelse(
  limma.results$logFC < 0 & limma.results$adj.P.Val <0.05, 'cyan', 
  ifelse(limma.results$logFC > 0  & limma.results$adj.P.Val <0.05, 'magenta', 'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'magenta'] <- 'upregulation'
names(keyvals)[keyvals == 'black'] <- 'NS'
names(keyvals)[keyvals == 'cyan'] <- 'downregulation'

EnhancedVolcano(limma.results,
                lab = rownames(limma.results),
                x = 'logFC',
                y = 'P.Value')
#Enhanced volcano plot
#to highlight some specific genes use the below one
specificgenes = c("COMT","3HAO","MAOM","AOFA","KAT3","STAT6","STA5B","FABP5","KMO","DGLB","FABPH"
,"TGFI1","SYWC","LHPP","1433S","LIPL","MGLL","I23O1", "PP2AA", "ACTH","SUMO2", "HOOK3")

p <- EnhancedVolcano(limma.results,
                     lab = rownames(limma.results),
                     x = 'logFC',
                     y = 'P.Value',
                     #selectLab = c(limma.results$Accession),
                     selectLab = specificgenes,
                     title= 'Human placental DEP (covariates controlled, without SANDY total50subj)',
                     #xlim = c(-4, 3),
                     #ylim = c(0, 5.5),
                     pCutoff= 0.05,
                     pointSize = 1,
                     labSize = 4.0,
                     #FCcutoff= 0.0,
                     colCustom = keyvals,
                     #shapeCustom = keyvals.shape,
                     #pointSize = c(ifelse(rownames(DEP.MS) %in% immunegene$Gene.Names & DEP.MS$Anova.p< 0.05, 4, 2)),
                     #shape = c(6, 6, 19, 16),
                     legendPosition = 'bottom',
                     borderColour = 'black',
                     legendLabSize = 10,
                     border = 'full',
                     shape = 'circle',
                     #boxedLabels = TRUE,
                     #labFace = 'bold',
                     #colGradient = c('yellow','blue', 'red'),
                     #borderColour = 'black',
                     #label_bquote(),
                     colAlpha = 4/5,
                     #scale_alpha_continuous(),
                     #borderWidth = 1.5,
                     #borderColour = 'black',
)

