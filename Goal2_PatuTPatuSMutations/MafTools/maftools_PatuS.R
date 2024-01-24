## ----eval=FALSE---------------------------------------------------------------
#  if (!require("BiocManager"))
#      install.packages("BiocManager")
#  BiocManager::install("maftools")

## ----loadlib, results='hide', message=FALSE-----------------------------------
library(maftools)

## ----readmaf------------------------------------------------------------------
#path to TCGA LAML MAF file
laml.maf = '/Directorio/PatuSOmicsSomaticMutationsMAFProfile.maf'

laml = read.maf(maf = laml.maf)

## ----mafobject----------------------------------------------------------------
#Typing laml shows basic summary of MAF file.
laml

## ----mafsummary, eval=FALSE---------------------------------------------------
#  #Shows sample summry.
getSampleSummary(laml)
#  #Shows gene summary.
getGeneSummary(laml)
#  #shows clinical data associated with samples
getClinicalData(laml)
#  #Shows all fields in MAF
getFields(laml)
#  #Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')

## ----summaryPlot,fig.height=4, fig.width=6------------------------------------
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

## ----titv, fig.height=3, fig.width=4.2, eval = T, fig.align='left'------------
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)

## ----drugInteractions, fig.height=3, fig.width=5------------------------------
dgi = drugInteractions(maf = laml, fontSize = 0.75)