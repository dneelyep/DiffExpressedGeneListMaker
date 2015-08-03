#####################################################
## THESE LINES NEED TO BE EXECUTED ONCE TO INSTALL ##
## THE CORRECT LIBRARIES -- THEY ARE COMMENTED OUT ##
## HERE SINCE THE LIBRARIES ARE ALREADY INSTALLED  ##
#####################################################

###########################################################
## DATA IS FROM Affymetrix Arabidopsis Gene 1.0 ST Array ##
###########################################################

##############################################
## INSTALL ALL OF THE BIOCONDUCTOR PACKAGES ##
##############################################
# TODO Convert all function definitions to use roxygen, and set up automatic
#      documentation generation, given that this is a library to be used by others.
installPackages <- function() {

  #   ##################################################################
  #   ## FUNCTION installPackages                                     ##
  #   ## Installs all packages required by this project on            ##
  #   ## this computer.                                               ##
  #   ##################################################################

  source("http://www.bioconductor.org/biocLite.R")
  
  print("Installing packages now...")
  
  biocLite("agcdf", dependencies=TRUE)
  biocLite("agprobe", dependencies=TRUE)
  biocLite("ath1121501.db", dependencies=TRUE)
  biocLite("ag.db", dependencies=TRUE)
  biocLite("pd.aragene.1.0.st", dependencies=TRUE)
  biocLite("affyQCReport", dependencies=TRUE)
  biocLite("arrayQualityMetrics", dependencies=TRUE)
  biocLite("hgu133plus2.db", dependencies=TRUE)
  biocLite("affyPLM", dependencies=TRUE)
  biocLite("RCytoscape", dependencies=TRUE)
  biocLite("categoryCompare", dependencies=TRUE)
  biocLite("AnnotationDbi", dependencies=TRUE)
  biocLite("GO.db", dependencies=TRUE)
  biocLite("limma", dependencies=TRUE)
  biocLite("puma")
  install.packages("ggplot2")
  install.packages("VennDiagram")
  install.packages("KEGG.db")
  
  print("...Done installing packages.")
}

installPackages()


# ###################################
# ## INCLUDE THE REQUIRED PACKAGES ##
# ###################################
# 
# require(org.At.eg.db)
# require(org.At.tair.db)
# require(ag.db)
# require(ath1121501.db)
# require(agcdf)
# require(agprobe)
# require(pd.aragene.1.0.st)
# require(oligo)
# require(limma)
# require(grid)
# require(VennDiagram)
# require(RCytoscape)
# require(graph)
# require(XMLRPC)
# require(categoryCompare)
# require(GO.db)
# require(KEGG.db)
# 
# #========================================================================================================
# getAnnotations <- function(dbPath, requestedID=c("ENSEMBL", "ENTREZ", "GENENAME", "RefSeq", "SYMBOL", "TAIR"), sepR="PrimeView") {
#   
#   ##################################################################
#   ## FUNCTION getAnnotations                                      ##
#   ## Anntoations are stored in Annotations/PrimeViewENSEMBL.annot ##
#   ## that have been created from primeview GEO GPL file           ##
#   ##                                                              ##
#   ## IN: path to the annotation database files                    ##
#   ##     list of annotations to add                               ##
#   ## OUT: Merged matrix of annotations containing all of the      ##
#   ##      specified annotations and the ID (affyProbeSetID)       ##
#   ##################################################################
#   
#   for(idType in requestedID) {
#     annFileName <- paste(dbPath, idType, sep = sepR);
#     annFileName <- paste(annFileName, ".annot", sep="");
#     
#     currAnnot <- read.table(annFileName, header=F, sep="\t", quote="")
#     names(currAnnot) <- c("ID", idType)  ## ID is Affy probeset name
#     if(!exists("allAnnot")) {
#       allAnnot<-currAnnot;               ## First pass through
#     }
#     else {
#       allAnnot<-merge(allAnnot, currAnnot, by="ID")
#     }
#   }
#   return(allAnnot)
# }
# #========================================================================================================
# 
# #========================================================================================================
# getComparisons <- function(fitData, comparison, adjMethod="BH", pValCut=1) {
#   #################################################################  
#   ## FUNCTION getComparisons                                     ##
#   ## DESCRIPTION: returns the comparison data from the fit data ##
#   ##              according to the desired comparison            ##
#   ## IN: fitData array (required)                                ##
#   ##             comparison to make (required)                   ##
#   ##             FDR adjustment method (BH by default)           ##
#   ##             PValue cutoff (1 by default)                    ##
#   ## OUT: table of comparison data, sorted by ID (affyProbeSet)  ##
#   #################################################################
#   allCompData <- topTable(fitData, comparison, number=Inf, adjust.method=adjMethod, p.value=pValCut)
#   allCompData <- allCompData[order(allCompData$ID), ]
#   return(allCompData)
# }
# #========================================================================================================
# 
# #========================================================================================================
# getDiffExpressedGenes <- function(allData, cutOffLabel="P.Value", cutOffLevel=0.05) {
#   
#   #################################################################  
#   ## FUNCTION getDiffExpressedGenes                              ##
#   ## DESCRIPTION:  Reads in a table of genes and their associated##
#   ## IN: allData (required)                                      ##
#   ##     cutOffLabel (P.Value by default) -- field label to make ##
#   ##     cutoffs for reporting genes                             ##
#   ##     cutOffLevel (0.05 by default) -- value by which genes   ##
#   ##     are included in the differentially expressed list       ##
#   ## OUT: de_all -- a list of differentially expressed genes,    ##
#   ##                sorted by their P value                      ##
#   #################################################################
#   
#   pLoc<-grep(cutOffLabel, names(allData), value = T)
#   de_all<-(allData[,pLoc[1]] <= cutOffLevel)
#   de_all<-allData[de_all, ]
#   de_all<-de_all[order(de_all$P.Value),]    #just sorts by P value -- has been cut by whatever label is present
#   return(de_all)
# }
# #========================================================================================================
# 
# #========================================================================================================
# findDEGenes_NoAnnot <- function (fData, comp, cutoffLabel="P.Value", cutoffLevel="0.05") {
#   #################################################################  
#   ## FUNCTION findDEGenes                                        ##
#   ## DESCRIPTION:                                                ##
#   ## IN:                                                         ##
#   ## OUT:                                                        ##
#   #################################################################
#   deData <- getComparisons(fitData=fData, comparison=comp)
#   deData <- getDiffExpressedGenes(allData=deData, cutOffLabel=cutoffLabel, cutOffLevel=cutoffLevel)
#   return(deData)
# }
# #========================================================================================================
# #========================================================================================================
# findALLGenes_NoAnnot <- function(fData, comp) {
#   #################################################################  
#   ## FUNCTION findALLGenes                                       ##
#   ## DESCRIPTION:                                                ##
#   ## IN:                                                         ##
#   ## OUT:                                                        ##
#   #################################################################
#   allGeneData <- getComparisons(fitData=fData, comparison=comp)
#   return(allGeneData)
# }
# #========================================================================================================
# #========================================================================================================
# findDEGenes <- function (fData, comp, cutoffLabel="P.Value", cutoffLevel="0.05", annotData) {
#   #################################################################  
#   ## FUNCTION findDEGenes                                        ##
#   ## DESCRIPTION:                                                ##
#   ## IN:                                                         ##
#   ## OUT:                                                        ##
#   #################################################################
#   deData <- getComparisons(fitData=fData, comparison=comp)
#   deData <- merge(deData, annotData, by="ID")
#   deData <- getDiffExpressedGenes(allData=deData, cutOffLabel=cutoffLabel, cutOffLevel=cutoffLevel)
#   return(deData)
# }
# #========================================================================================================
# 
# #========================================================================================================
# findALLGenes <- function(fData, comp, annotData) {
#   #################################################################  
#   ## FUNCTION findALLGenes                                       ##
#   ## DESCRIPTION:                                                ##
#   ## IN:                                                         ##
#   ## OUT:                                                        ##
#   #################################################################
#   allGeneData <- getComparisons(fitData=fData, comparison=comp)
#   allGeneData <- merge(allGeneData, annotData, by="ID")
#   return(allGeneData)
# }
# #========================================================================================================
# 
# #========================================================================================================
# findDEGs <- function(fit.b, comp1, comp2, cutoffLabel, cutoffLevel, annotData) {
#   #################################################################  
#   ## FUNCTION findDEGs                                           ##
#   ## DESCRIPTION:                                                ##
#   ## IN:                                                         ##
#   ## OUT:                                                        ##
#   #################################################################
#   tmp = paste(comp1, "-", comp2)
#   DEGs = findDEGenes(fData=fit.b, comp=tmp, cutoffLabel=cutoffLabel, cutoffLevel=cutoffLevel, annotData=annotData)
#   return(DEGs)
# }
# #========================================================================================================
# 
# #========================================================================================================
# findALL <- function(fit.b, comp1, comp2, allAnn) {
#   #################################################################  
#   ## FUNCTION findALL                                            ##
#   ## DESCRIPTION:                                                ##
#   ## IN:                                                         ##
#   ## OUT:                                                        ##
#   #################################################################
#   tmp = paste(comp1, "-", comp2)
#   ALLGenes = findALLGenes(fData=fit.b, comp=tmp, annotData=allAnn)
#   return(ALLGenes)
# }
# #========================================================================================================
# 
# #========================================================================================================
# createFile<-function(GeneList, baseType, comp1, comp2) {
#   #################################################################  
#   ## FUNCTION createFile                                         ##
#   ## DESCRIPTION:                                                ##
#   ## IN:                                                         ##
#   ## OUT:                                                        ##
#   #################################################################
#   outFileName = paste(baseType, "_", comp1, "_vs_", comp2, ".tab", sep="")
#   write.table(GeneList, file=outFileName, sep="\t", col.names=T, row.names=F)
# }
# #========================================================================================================
# 
# #========================================================================================================
# constructThreeVennSet <- function(geneSet1, geneSet2, geneSet3, geneSetLabels, fillColors, fName) {
#   
#   ######################################################
#   ## Constructs a venn diagram containing information ##
#   ## on the number of items in three sets             ##
#   ######################################################
#   
#   #####################################################
#   ## Remove duplicate entries and NAs from gene sets ##
#   #####################################################
#   geneSet1 <- unique(geneSet1)
#   geneSet1 <- geneSet1[which(geneSet1 != "NA")]
#   geneSet2 <- unique(geneSet2)
#   geneSet2 <- geneSet2[which(geneSet2 != "NA")]
#   geneSet3 <- unique(geneSet3)
#   geneSet3 <- geneSet3[which(geneSet3 != "NA")]
#   
#   Set12 <- intersect(geneSet1, geneSet2)
#   Set23 <- intersect(geneSet2, geneSet3)
#   Set13 <- intersect(geneSet1, geneSet3)
#   Set123 <- intersect(Set12, geneSet3)
#   
#   Venn.Root <- draw.triple.venn(length(geneSet1), length(geneSet2), length(geneSet3), length(Set12), length(Set23), 
#                                 length(Set13), length(Set123), category=geneSetLabels,
#                                 fill=fillColors)
#   tiff(filename=fName, compression="lzw")
#   grid.draw(Venn.Root)  
#   dev.off()
#   return(Set123)
# }
# #========================================================================================================
# 
# #========================================================================================================
# getCCAnnotationTableResults <- function(annotTable, annotationDescription){
#   annotTable <- merge(annotTable, annotationDescription, by="ID")
#   annotTable <- annotTable[order(annotTable$Pvalue),]
#   return(annotTable)
# }
# #========================================================================================================
# 
# #========================================================================================================
# runQAQC <- function(celData, genePS, geneCore, rmaDat_genePS, rmaDat_geneCore) {
#   ############################################
#   ## View raw image and non-normalized data ##
#   ############################################
#   image(celData, 1)
#   MAplot(celData[,1:4], pairs=TRUE)
#   boxplot(celData)
#   hist(celData)
#   boxplot(genePS)
#   boxplot(geneCore)
#   hist(genePS)
#   hist(geneCore)
#   boxplot(rmaDat_genePS)
#   boxplot(rmaDat_geneCore)
#   
#   ##############################################
#   ## LOOK AT CORRELATIONS BETWEEN EXPERIMENTS ##
#   ##############################################
#   corDat_rmaDat_genePS   <- cor(rmaDat_genePS)
#   corDat_rmaDat_geneCore <- cor(rmaDat_geneCore)
#   apply(corDat_rmaDat_genePS, 2, min)
#   apply(corDat_rmaDat_geneCore, 2, min)
# }
# #========================================================================================================
# 
# #========================================================================================================
# getDiffExpressedGenesFCcutoff <- function(allData, cutoffLabel="logFC", cutoffLevel=2.0) {
#   
#   #################################################################  
#   ## FUNCTION getDiffExpressedGenes                              ##
#   ## DESCRIPTION:  Reads in a table of genes and their associated##
#   ## IN: allData (required)                                      ##
#   ##     cutOffLabel (P.Value by default) -- field label to make ##
#   ##     cutoffs for reporting genes                             ##
#   ##     cutOffLevel (0.05 by default) -- value by which genes   ##
#   ##     are included in the differentially expressed list       ##
#   ## OUT: de_all -- a list of differentially expressed genes,    ##
#   ##                sorted by their P value                      ##
#   #################################################################
#   
#   pLoc<-grep(cutoffLabel, names(allData), value = T)
#   de_all<-(abs(allData[,pLoc[1]]) >= cutoffLevel)
#   de_all<-allData[de_all, ]
#   de_all<-de_all[order(de_all$P.Value),]    #just sorts by P value -- has been cut by whatever label is present
#   return(de_all)
# }
# #========================================================================================================
# 
# 
# 
# #========================================================================================================
# ######################################### MAIN PROGRAM ##################################################
# #========================================================================================================
# 
# 
# #############################
# # Read in the cel files for #
# # analysis using the affy   #
# # routines                  #
# #############################
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/")
# dir(datadir)
# setwd(datadir)
# celfiles<-list.celfiles()
# celData <- read.celfiles(files=celfiles)
# 
# allAnn_probe <- getAnnotations(dbPath="E:/Arabidopsis_probe_annotations/")
# allAnn_core  <- getAnnotations(dbPath="E:/Arabidopsis_core_annotations/")
# 
# ###################################################
# ## Retrieve the sample names and remove the      ##
# ## prefix and suffix, leaving only the treatment ##
# ## and tissue types                              ##
# ###################################################
# sNames <- sampleNames(celData)
# sNames <- gsub('^SS_', '', sNames)
# sNames <- gsub('_1.CEL', '#1', sNames)
# sNames <- gsub('_2.CEL', '#2', sNames)
# sNames <- gsub('_3.CEL', '#3', sNames)
# sampleNames(celData) <- sNames
# sampleNames(celData)
# 
# ###################################
# ## Now create a phenodata object ##
# ###################################
# pd <- do.call(rbind, strsplit(sNames, '#'))
# pd <- as.data.frame(pd)
# names(pd) <- c('treament_and_tissue', 'replicate')
# pd <- new('AnnotatedDataFrame', data=pd)
# sampleNames(pd) <- sNames
# phenoData(celData) <- pd
# rm(pd, sNames)
# 
# ############################
# ## Probe level model fits ##
# ############################
# plmFit_PS   <- fitProbeLevelModel(celData, target="probeset")
# plmFit_CORE <- fitProbeLevelModel(celData, target="core")
# 
# genePS   <- rma(celData, target="probeset")  ## Probeset level
# geneCore <- rma(celData, target="core")    ## Gene level
# 
# featureData(genePS)   <- getNetAffx(genePS, "probeset")
# featureData(geneCore) <- getNetAffx(geneCore, "transcript")
# 
# targets <- readTargets("phenoData.txt", path=datadir, sep="\t", row.names="filename")
# rmaDat_genePS   <- exprs(genePS)
# rmaDat_geneCore <- exprs(geneCore)
# 
# #########################################
# ## NOW WRITE THE DATA MATRICES FOR GEO ##
# #########################################
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/Results/GEO/")
# dir(datadir)
# setwd(datadir)
# 
# write.table(rmaDat_genePS,   file="GEO-ProbeLevelExpressionGPL18349.txt", sep="\t", col.names=T, row.names=T)
# write.table(rmaDat_geneCore, file="GEO-TranscriptLevelExpressionGPL17416.txt", sep="\t", col.names=T, row.names=T)
# 
# 
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/")
# dir(datadir)
# setwd(datadir)
# # runQAQC (celData=celData, genePS=genePS, geneCore=geneCore, rmaDat_genePS=rmaDat_genePS, rmaDat_geneCore=rmaDat_geneCore)
# 
# ############################################
# ## DESIGN AT THE PROBE SET LEVEL ANALYSIS ##
# ############################################
# f<-factor(targets$tissue_and_treatment, levels=unique(targets$tissue_and_treatment))
# design<-model.matrix(~0 + f)
# colnames(design)<-levels(f)
# fit<-lmFit(exprs(genePS), design)
# 
# #######################################
# ## DESIGN AT THE GENE LEVEL ANALYSIS ##
# #######################################
# f_core<-factor(targets$tissue_and_treatment, levels=unique(targets$tissue_and_treatment))
# design_core<-model.matrix(~0 + f_core)
# colnames(design_core)<-levels(f_core)
# fit_core<-lmFit(exprs(geneCore), design_core)
# unique(targets$tissue_and_treatment)
# phenoData(celData)
# 
# ROOT_contrast.matrix <- makeContrasts(TISSUE_AND_TREATMENT_Root_20mM       - TISSUE_AND_TREATMENT_Root_1.25mM,
#                                       TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_1.25mM,
#                                       TISSUE_AND_TREATMENT_Root_P          - TISSUE_AND_TREATMENT_Root_1.25mM,
#                                       TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_20mM,
#                                       TISSUE_AND_TREATMENT_Root_P          - TISSUE_AND_TREATMENT_Root_20mM,
#                                       TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_P,
#                                       levels=design_core)
# 
# SHOOT_contrast.matrix <- makeContrasts(TISSUE_AND_TREATMENT_Shoot_20mM       - TISSUE_AND_TREATMENT_Shoot_1.25mM,
#                                        TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_1.25mM,
#                                        TISSUE_AND_TREATMENT_Shoot_P          - TISSUE_AND_TREATMENT_Shoot_1.25mM,
#                                        TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_20mM,
#                                        TISSUE_AND_TREATMENT_Shoot_P          - TISSUE_AND_TREATMENT_Shoot_20mM,
#                                        TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_P,
#                                        levels=design_core)
# 
# ROOT_SHOOT_contrast.matrix <- makeContrasts(TISSUE_AND_TREATMENT_Shoot_20mM - TISSUE_AND_TREATMENT_Root_20mM,
#                                             TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Root_20mMS6Pct,
#                                             TISSUE_AND_TREATMENT_Shoot_1.25mM - TISSUE_AND_TREATMENT_Root_1.25mM,
#                                             TISSUE_AND_TREATMENT_Shoot_P - TISSUE_AND_TREATMENT_Root_P,
#                                             levels=design_core)
# 
# 
# #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
# #----------------------------------------------------------------------------------------------------------------------------------#
# #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv#
# 
# ################################
# ## ROOT DIFFERENTIAL ANALYSIS ##
# ################################
# 
# fit_ROOT   <-contrasts.fit(fit_core, ROOT_contrast.matrix)
# fit_ROOT.b <- eBayes(fit_ROOT)
# de_TREATED_ROOT_20mM_vs_ROOT_1.25mM     <- findDEGenes (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mM - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core)
# de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_ALL <- findALLGenes(fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mM - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core)
# de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_adj <- findDEGenes (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mM - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_adj_FC <- getDiffExpressedGenesFCcutoff (de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM     <- findDEGenes (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core)
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_ALL <- findALLGenes(fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core)
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_adj <- findDEGenes (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_adj_FC <- getDiffExpressedGenesFCcutoff (de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_ROOT_P_vs_ROOT_1.25mM     <- findDEGenes  (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_P - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core)
# de_TREATED_ROOT_P_vs_ROOT_1.25mM_ALL <- findALLGenes (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_P - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core)
# de_TREATED_ROOT_P_vs_ROOT_1.25mM_adj <- findDEGenes  (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_P - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_ROOT_P_vs_ROOT_1.25mM_adj_FC <- getDiffExpressedGenesFCcutoff (de_TREATED_ROOT_P_vs_ROOT_1.25mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_20mM     <- findDEGenes  (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_20mM", annotData=allAnn_core)
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_20mM_ALL <- findALLGenes (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_20mM", annotData=allAnn_core)
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_20mM_adj <- findDEGenes  (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_20mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_20mM_adj_FC <- getDiffExpressedGenesFCcutoff (de_TREATED_ROOT_20mMS6Pct_vs_ROOT_20mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_ROOT_P_vs_ROOT_20mM     <- findDEGenes  (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_P - TISSUE_AND_TREATMENT_Root_20mM", annotData=allAnn_core)
# de_TREATED_ROOT_P_vs_ROOT_20mM_ALL <- findALLGenes (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_P - TISSUE_AND_TREATMENT_Root_20mM", annotData=allAnn_core)
# de_TREATED_ROOT_P_vs_ROOT_20mM_adj <- findDEGenes  (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_P - TISSUE_AND_TREATMENT_Root_20mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_ROOT_P_vs_ROOT_20mM_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_ROOT_P_vs_ROOT_20mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_P     <- findDEGenes  (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_P", annotData=allAnn_core)
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_P_ALL <- findALLGenes (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_P", annotData=allAnn_core)
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_P_adj <- findDEGenes  (fData=fit_ROOT.b, comp="TISSUE_AND_TREATMENT_Root_20mMS6Pct - TISSUE_AND_TREATMENT_Root_P", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_ROOT_20mMS6Pct_vs_ROOT_P_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_P_adj, cutoffLevel=2.0)
# 
# ##############################################
# ## Write out differentially expressed genes ##
# ##############################################
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/Results/Differentially Expressed Genes/ROOT/")
# dir(datadir)
# setwd(datadir)
# write.table(de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_adj,     file="ROOT_20mMvs1.25mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_adj,file="ROOT_20mMS6Pctvs1.25mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_P_vs_ROOT_1.25mM_adj,        file="ROOT_Pvs1.25mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_20mM_adj,  file="ROOT_20mMS6Pctvs20mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_P_vs_ROOT_20mM_adj,          file="ROOT_Pvs20mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_P_adj,     file="ROOT_20mMS6PctvsP_adj.txt", sep="\t", col.names=T, row.names=F)
# 
# RootAllThree <- constructThreeVennSet(geneSet1=de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_adj$ENTREZ,
#                                       geneSet2=de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_adj$ENTREZ,
#                                       geneSet3=de_TREATED_ROOT_P_vs_ROOT_1.25mM_adj$ENTREZ,
#                                       geneSetLabels=c("20mM", "20mM-S6%", "P-"), fillColors=c("red", "green", "blue"),
#                                       fName="venn_root_pVal.tiff")
# 
# write.table(de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_adj_FC,     file="ROOT_20mMvs1.25mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_adj_FC,file="ROOT_20mMS6Pctvs1.25mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_P_vs_ROOT_1.25mM_adj_FC,        file="ROOT_Pvs1.25mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_20mM_adj_FC,  file="ROOT_20mMS6Pctvs20mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_P_vs_ROOT_20mM_adj_FC,          file="ROOT_Pvs20mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_P_adj_FC,     file="ROOT_20mMS6PctvsP_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# 
# constructThreeVennSet(geneSet1=de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_adj_FC$ENTREZ,
#                       geneSet2=de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_adj_FC$ENTREZ,
#                       geneSet3=de_TREATED_ROOT_P_vs_ROOT_1.25mM_adj_FC$ENTREZ,
#                       geneSetLabels=c("20mM", "20mM-S6%", "P-"), fillColors=c("red", "green", "blue"),
#                       fName="venn_root_pVal_FC.tiff")
# 
# 
# #########################
# ## Write out all genes ##
# #########################
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/Results/All Genes/ROOT/")
# dir(datadir)
# setwd(datadir)
# write.table(de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_ALL,     file="ROOT_20mMvs1.25mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_ALL,file="ROOT_20mMS6Pctvs1.25mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_P_vs_ROOT_1.25mM_ALL,        file="ROOT_Pvs1.25mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_20mM_ALL,  file="ROOT_20mMS6Pctvs20mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_P_vs_ROOT_20mM_ALL,          file="ROOT_Pvs20mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_P_ALL,     file="ROOT_20mMS6PctvsP_ALL.txt", sep="\t", col.names=T, row.names=F)
# 
# 
# 
# 
# ###########################################
# ## NOW DO categoryCompare for Enrichment ##
# ###########################################
# gUniverse <- unique(de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_ALL$TAIR)
# 
# table_ROOT_20mM_vs_ROOT_1.25mM <- unlist(de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_adj$TAIR)
# list_ROOT_20mM_vs_ROOT_1.25mM  <- list(genes=table_ROOT_20mM_vs_ROOT_1.25mM, universe=gUniverse, annotation='org.At.tair.db')
# 
# table_ROOT_20mMS6Pct_vs_ROOT_1.25mM <- unlist(de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_adj$TAIR)
# list_ROOT_20mMS6Pct_vs_ROOT_1.25mM  <- list(genes=table_ROOT_20mMS6Pct_vs_ROOT_1.25mM, universe=gUniverse, annotation='org.At.tair.db')
# 
# table_ROOT_P_vs_ROOT_1.25mM <- unlist(de_TREATED_ROOT_P_vs_ROOT_1.25mM_adj$TAIR)
# list_ROOT_P_vs_ROOT_1.25mM  <- list(genes=table_ROOT_P_vs_ROOT_1.25mM, universe=gUniverse, annotation='org.At.tair.db')
# 
# geneLists_ROOT <- list(T_20mM=list_ROOT_20mM_vs_ROOT_1.25mM, T_20mMS6Pct=list_ROOT_20mMS6Pct_vs_ROOT_1.25mM, T_P=list_ROOT_P_vs_ROOT_1.25mM)
# geneLists_ROOT <- new('ccGeneList', geneLists_ROOT, ccType=c('BP', 'KEGG'))
# 
# geneLists_ROOT11 <- list(T_20mM=list_ROOT_20mM_vs_ROOT_1.25mM, T_20mMS6Pct=list_ROOT_20mMS6Pct_vs_ROOT_1.25mM, T_P=list_ROOT_P_vs_ROOT_1.25mM)
# geneLists_ROOT11 <- new('ccGeneList', geneLists_ROOT11, ccType=c('BP', 'KEGG'))
# 
# fdr(geneLists_ROOT) <- 0
# enrichLists_ROOT <- ccEnrich(geneLists_ROOT)
# pvalueCutoff(enrichLists_ROOT$BP) <- 0.001
# 
# ccOpts_ROOT <- new('ccOptions', listNames = names(geneLists_ROOT), outType='none')
# ccResults_ROOT <- ccCompare(enrichLists_ROOT, ccOpts_ROOT)
# #cw_ROOT.BP <- ccOutCyt(ccResults_ROOT$BP, ccOpts_ROOT)
# #breakEdges(cw_ROOT.BP, 0.8)
# 
# #########################################################
# ## Retrieve Description information for annotation IDs ##
# #########################################################
# goDescAnnot <- read.table("E:/GO Annotations/GOtoDesc.dat", header=F, sep="\t", quote="")
# names(goDescAnnot) <- c("ID", "GO Term Description")
# 
# KEGGDescAnnot <- read.table("E:/KEGG Annotations/KEGGtoDesc.dat", header=F, sep="\t", quote="", colClasses=c(rep("factor", 2)))
# names(KEGGDescAnnot) <- c("ID", "KEGG Term Description")
# 
# ##################################################
# ## Gene Ontology Biological Process Annotations ##
# ##################################################
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/Results/Gene Ontology/ROOT/")
# dir(datadir)
# setwd(datadir)
# 
# GO_BP_ROOT_20mM      <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT[[1]][[1]]), annotationDescription=goDescAnnot)
# GO_BP_ROOT_20mMS6Pct <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT[[1]][[2]]), annotationDescription=goDescAnnot)
# GO_BP_ROOT_P         <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT[[1]][[3]]), annotationDescription=goDescAnnot)
# 
# write.table(GO_BP_ROOT_20mM,      file="GO-BP-ROOT_20mM.txt", sep="\t", col.names=T, row.names=F)
# write.table(GO_BP_ROOT_20mMS6Pct, file="GO-BP-ROOT_20mMS6Pct.txt", sep="\t", col.names=T, row.names=F)
# write.table(GO_BP_ROOT_P,         file="GO-BP-ROOT_P.txt", sep="\t", col.names=T, row.names=F)
# 
# ######################
# ## KEGG ANNOTATIONS ##
# ######################
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/Results/KEGG Pathways/ROOT/")
# dir(datadir)
# setwd(datadir)
# KEGG_ROOT_20mM      <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT[[2]][[1]]), annotationDescription=KEGGDescAnnot)
# KEGG_ROOT_20mMS6Pct <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT[[2]][[2]]), annotationDescription=KEGGDescAnnot)
# KEGG_ROOT_P         <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT[[2]][[3]]), annotationDescription=KEGGDescAnnot)
# 
# write.table(KEGG_ROOT_20mM,      file="KEGG-ROOT_20mM.txt", sep="\t", col.names=T, row.names=F)
# write.table(KEGG_ROOT_20mMS6Pct, file="KEGG-ROOT_20mMS6Pct.txt", sep="\t", col.names=T, row.names=F)
# write.table(KEGG_ROOT_P,         file="KEGG-ROOT_P.txt", sep="\t", col.names=T, row.names=F)
# 
# 
# #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
# #----------------------------------------------------------------------------------------------------------------------------------#
# #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv#
# 
# #################################
# ## SHOOT DIFFERENTIAL ANALYSIS ##
# #################################
# fit_SHOOT   <- contrasts.fit(fit_core, SHOOT_contrast.matrix)
# fit_SHOOT.b <- eBayes(fit_SHOOT)
# de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM     <- findDEGenes (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mM - TISSUE_AND_TREATMENT_Shoot_1.25mM", annotData=allAnn_core)
# de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_ALL <- findALLGenes(fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mM - TISSUE_AND_TREATMENT_Shoot_1.25mM", annotData=allAnn_core)
# de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_adj <- findDEGenes (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mM - TISSUE_AND_TREATMENT_Shoot_1.25mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_adj_FC <- getDiffExpressedGenesFCcutoff (de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM     <- findDEGenes (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_1.25mM", annotData=allAnn_core)
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_ALL <- findALLGenes(fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_1.25mM", annotData=allAnn_core)
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_adj <- findDEGenes (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_1.25mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_SHOOT_P_vs_SHOOT_1.25mM     <- findDEGenes  (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_P - TISSUE_AND_TREATMENT_Shoot_1.25mM", annotData=allAnn_core)
# de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_ALL <- findALLGenes (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_P - TISSUE_AND_TREATMENT_Shoot_1.25mM", annotData=allAnn_core)
# de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_adj <- findDEGenes  (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_P - TISSUE_AND_TREATMENT_Shoot_1.25mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_20mM     <- findDEGenes  (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_20mM", annotData=allAnn_core)
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_20mM_ALL <- findALLGenes (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_20mM", annotData=allAnn_core)
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_20mM_adj <- findDEGenes  (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_20mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_20mM_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_20mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_SHOOT_P_vs_SHOOT_20mM     <- findDEGenes  (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_P - TISSUE_AND_TREATMENT_Shoot_20mM", annotData=allAnn_core)
# de_TREATED_SHOOT_P_vs_SHOOT_20mM_ALL <- findALLGenes (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_P - TISSUE_AND_TREATMENT_Shoot_20mM", annotData=allAnn_core)
# de_TREATED_SHOOT_P_vs_SHOOT_20mM_adj <- findDEGenes  (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_P - TISSUE_AND_TREATMENT_Shoot_20mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_SHOOT_P_vs_SHOOT_20mM_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_SHOOT_P_vs_SHOOT_20mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_P     <- findDEGenes  (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_P", annotData=allAnn_core)
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_P_ALL <- findALLGenes (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_P", annotData=allAnn_core)
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_P_adj <- findDEGenes  (fData=fit_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Shoot_P", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_P_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_P_adj, cutoffLevel=2.0)
# 
# 
# ##############################################
# ## Write out differentially expressed genes ##
# ##############################################
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/Results/Differentially Expressed Genes/SHOOT/")
# dir(datadir)
# setwd(datadir)
# write.table(de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_adj,file="SHOOT_20mMvs1.25mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_adj,file="SHOOT_20mMS6Pctvs1.25mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_adj,file="SHOOT_Pvs1.25mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_20mM_adj,file="SHOOT_20mMS6Pctvs20mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_P_vs_SHOOT_20mM_adj,file="SHOOT_Pvs20mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_P_adj,file="SHOOT_20mMS6PctvsP_adj.txt", sep="\t", col.names=T, row.names=F)
# 
# write.table(de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_adj_FC,file="SHOOT_20mMvs1.25mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_adj_FC,file="SHOOT_20mMS6Pctvs1.25mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_adj_FC,file="SHOOT_Pvs1.25mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_20mM_adj_FC,file="SHOOT_20mMS6Pctvs20mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_P_vs_SHOOT_20mM_adj_FC,file="SHOOT_Pvs20mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_P_adj_FC,file="SHOOT_20mMS6PctvsP_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# 
# ############################
# ## Construct VENN Diagram ##
# ############################
# ShootAllThree <- constructThreeVennSet(geneSet1=de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_adj$ENTREZ,
#                                        geneSet2=de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_adj$ENTREZ,
#                                        geneSet3=de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_adj$ENTREZ,
#                                        geneSetLabels=c("20mM", "20mM-S6%", "P-"), fillColors=c("red", "green", "blue"),
#                                        fName="venn_SHOOT_pval.tiff")
# 
# constructThreeVennSet(geneSet1=de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_adj_FC$ENTREZ,
#                       geneSet2=de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_adj_FC$ENTREZ,
#                       geneSet3=de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_adj_FC$ENTREZ,
#                       geneSetLabels=c("20mM", "20mM-S6%", "P-"), fillColors=c("red", "green", "blue"),
#                       fName="venn_SHOOT_pval_FC.tiff")
# 
# 
# #########################
# ## Write out all genes ##
# #########################
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/Results/All Genes/SHOOT/")
# dir(datadir)
# setwd(datadir)
# write.table(de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_ALL,file="SHOOT_20mMvs1.25mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_ALL,file="SHOOT_20mMS6Pctvs1.25mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_ALL,file="SHOOT_Pvs1.25mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_20mM_ALL,file="SHOOT_20mMS6Pctvs20mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_P_vs_SHOOT_20mM_ALL,file="SHOOT_Pvs20mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_P_ALL,file="SHOOT_20mMS6PctvsP_ALL.txt", sep="\t", col.names=T, row.names=F)
# 
# 
# 
# 
# ###########################################
# ## NOW DO categoryCompare for Enrichment ##
# ###########################################
# table_SHOOT_20mM_vs_SHOOT_1.25mM <- unlist(de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_adj$TAIR)
# list_SHOOT_20mM_vs_SHOOT_1.25mM  <- list(genes=table_SHOOT_20mM_vs_SHOOT_1.25mM, universe=gUniverse, annotation='org.At.tair.db')
# 
# table_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM <- unlist(de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_adj$TAIR)
# list_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM  <- list(genes=table_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM, universe=gUniverse, annotation='org.At.tair.db')
# 
# table_SHOOT_P_vs_SHOOT_1.25mM <- unlist(de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_adj$TAIR)
# list_SHOOT_P_vs_SHOOT_1.25mM  <- list(genes=table_SHOOT_P_vs_SHOOT_1.25mM, universe=gUniverse, annotation='org.At.tair.db')
# 
# geneLists_SHOOT <- list(T_20mM=list_SHOOT_20mM_vs_SHOOT_1.25mM, T_20mMS6Pct=list_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM, T_P=list_SHOOT_P_vs_SHOOT_1.25mM)
# geneLists_SHOOT <- new('ccGeneList', geneLists_SHOOT, ccType=c('BP', 'KEGG'))
# 
# fdr(geneLists_SHOOT) <- 0
# enrichLists_SHOOT <- ccEnrich(geneLists_SHOOT)
# pvalueCutoff(enrichLists_SHOOT$BP) <- 0.001
# 
# ccOpts_SHOOT <- new('ccOptions', listNames = names(geneLists_SHOOT), outType='none')
# ccResults_SHOOT <- ccCompare(enrichLists_SHOOT, ccOpts_SHOOT)
# cw_SHOOT.BP <- ccOutCyt(ccResults_SHOOT$BP, ccOpts_SHOOT)
# breakEdges(cw_SHOOT.BP, 0.8)
# 
# ##################################################
# ## Gene Ontology Biological Process Annotations ##
# ##################################################
# 
# GO_BP_SHOOT_20mM      <- getCCAnnotationTableResults(annotTable=summary(enrichLists_SHOOT[[1]][[1]]), annotationDescription=goDescAnnot)
# GO_BP_SHOOT_20mMS6Pct <- getCCAnnotationTableResults(annotTable=summary(enrichLists_SHOOT[[1]][[2]]), annotationDescription=goDescAnnot)
# GO_BP_SHOOT_P         <- getCCAnnotationTableResults(annotTable=summary(enrichLists_SHOOT[[1]][[3]]), annotationDescription=goDescAnnot)
# 
# write.table(GO_BP_SHOOT_20mM,      file="GO-BP-SHOOT_20mM.txt", sep="\t", col.names=T, row.names=F)
# write.table(GO_BP_SHOOT_20mMS6Pct, file="GO-BP-SHOOT_20mMS6Pct.txt", sep="\t", col.names=T, row.names=F)
# write.table(GO_BP_SHOOT_P,         file="GO-BP-SHOOT_P.txt", sep="\t", col.names=T, row.names=F)
# 
# ######################
# ## KEGG ANNOTATIONS ##
# ######################
# KEGG_SHOOT_20mM      <- getCCAnnotationTableResults(annotTable=summary(enrichLists_SHOOT[[2]][[1]]), annotationDescription=KEGGDescAnnot)
# KEGG_SHOOT_20mMS6Pct <- getCCAnnotationTableResults(annotTable=summary(enrichLists_SHOOT[[2]][[2]]), annotationDescription=KEGGDescAnnot)
# KEGG_SHOOT_P         <- getCCAnnotationTableResults(annotTable=summary(enrichLists_SHOOT[[2]][[3]]), annotationDescription=KEGGDescAnnot)
# 
# write.table(KEGG_SHOOT_20mM,      file="KEGG-SHOOT_20mM.txt", sep="\t", col.names=T, row.names=F)
# write.table(KEGG_SHOOT_20mMS6Pct, file="KEGG-SHOOT_20mMS6Pct.txt", sep="\t", col.names=T, row.names=F)
# write.table(KEGG_SHOOT_P,         file="KEGG-SHOOT_P.txt", sep="\t", col.names=T, row.names=F)
# 
# 
# 
# 
# ### WRITE OUT THE ENTREZ IDs FOR ALL GENES
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/Results/Heatmap/")
# dir(datadir)
# setwd(datadir)
# 
# Root_PVal_FC_UNION <- union(de_TREATED_ROOT_20mM_vs_ROOT_1.25mM_adj_FC$ENTREZ,
#                             de_TREATED_ROOT_20mMS6Pct_vs_ROOT_1.25mM_adj_FC$ENTREZ)
# Root_PVal_FC_UNION <- union(Root_PVal_FC_UNION, de_TREATED_ROOT_P_vs_ROOT_1.25mM_adj_FC$ENTREZ)
# write.table(Root_PVal_FC_UNION,     file="ROOT_DEGs_PVal_FC.txt", sep="\t", col.names=F, row.names=F)
# Shoot_PVal_FC_UNION <- union(de_TREATED_SHOOT_20mM_vs_SHOOT_1.25mM_adj_FC$ENTREZ,
#                              de_TREATED_SHOOT_20mMS6Pct_vs_SHOOT_1.25mM_adj_FC$ENTREZ)
# Shoot_PVal_FC_UNION <- union(Shoot_PVal_FC_UNION, de_TREATED_SHOOT_P_vs_SHOOT_1.25mM_adj_FC$ENTREZ)
# write.table(Shoot_PVal_FC_UNION,     file="SHOOT_DEGs_PVal_FC.txt", sep="\t", col.names=F, row.names=F)
# 
# #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
# #----------------------------------------------------------------------------------------------------------------------------------#
# #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv#
# 
# ##########################################
# ## ROOT AND SHOOT DIFFERENTIAL ANALYSIS ##
# ##########################################
# 
# 
# fit_ROOT_SHOOT   <- contrasts.fit(fit_core, ROOT_SHOOT_contrast.matrix)
# fit_ROOT_SHOOT.b <- eBayes(fit_ROOT_SHOOT)
# 
# de_TREATED_ROOT_SHOOT_20mM     <- findDEGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mM - TISSUE_AND_TREATMENT_Root_20mM", annotData=allAnn_core)
# de_TREATED_ROOT_SHOOT_20mM_ALL <- findALLGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mM - TISSUE_AND_TREATMENT_Root_20mM", annotData=allAnn_core)
# de_TREATED_ROOT_SHOOT_20mM_adj <- findDEGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mM - TISSUE_AND_TREATMENT_Root_20mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_ROOT_SHOOT_20mM_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_ROOT_SHOOT_20mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_ROOT_SHOOT_20mMS6Pct     <- findDEGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Root_20mMS6Pct", annotData=allAnn_core)
# de_TREATED_ROOT_SHOOT_20mMS6Pct_ALL <- findALLGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Root_20mMS6Pct", annotData=allAnn_core)
# de_TREATED_ROOT_SHOOT_20mMS6Pct_adj <- findDEGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_20mMS6Pct - TISSUE_AND_TREATMENT_Root_20mMS6Pct", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_ROOT_SHOOT_20mMS6Pct_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_ROOT_SHOOT_20mMS6Pct_adj, cutoffLevel=2.0)
# 
# de_TREATED_ROOT_SHOOT_1.25mM     <- findDEGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_1.25mM - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core)
# de_TREATED_ROOT_SHOOT_1.25mM_ALL <- findALLGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_1.25mM - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core)
# de_TREATED_ROOT_SHOOT_1.25mM_adj <- findDEGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_1.25mM - TISSUE_AND_TREATMENT_Root_1.25mM", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_ROOT_SHOOT_1.25mM_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_ROOT_SHOOT_1.25mM_adj, cutoffLevel=2.0)
# 
# de_TREATED_ROOT_SHOOT_P     <- findDEGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_P - TISSUE_AND_TREATMENT_Root_P", annotData=allAnn_core)
# de_TREATED_ROOT_SHOOT_P_ALL <- findALLGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_P - TISSUE_AND_TREATMENT_Root_P", annotData=allAnn_core)
# de_TREATED_ROOT_SHOOT_P_adj <- findDEGenes(fData=fit_ROOT_SHOOT.b, comp="TISSUE_AND_TREATMENT_Shoot_P - TISSUE_AND_TREATMENT_Root_P", annotData=allAnn_core, cutoffLabel="adj.P.Val", cutoffLevel=0.05)
# de_TREATED_ROOT_SHOOT_P_adj_FC <- getDiffExpressedGenesFCcutoff(de_TREATED_ROOT_SHOOT_P_adj, cutoffLevel=2.0)
# 
# 
# ##############################################
# ## Write out differentially expressed genes ##
# ##############################################
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/Results/Differentially Expressed Genes/ROOT-SHOOT/")
# dir(datadir)
# setwd(datadir)
# write.table(de_TREATED_ROOT_SHOOT_20mM_adj,file="ROOT_SHOOT_20mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_SHOOT_20mMS6Pct_adj,file="ROOT_SHOOT_20mMS6Pct_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_SHOOT_1.25mM_adj,file="ROOT_SHOOT_1.25mM_adj.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_SHOOT_P_adj,file="ROOT_SHOOT_P_adj.txt", sep="\t", col.names=T, row.names=F)
# 
# write.table(de_TREATED_ROOT_SHOOT_20mM_adj_FC,file="ROOT_SHOOT_20mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_SHOOT_20mMS6Pct_adj_FC,file="ROOT_SHOOT_20mMS6Pct_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_SHOOT_1.25mM_adj_FC,file="ROOT_SHOOT_1.25mM_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_SHOOT_P_adj_FC,file="ROOT_SHOOT_P_adj_FC.txt", sep="\t", col.names=T, row.names=F)
# 
# #########################
# ## Write out all genes ##
# #########################
# datadir <- file.path("E:/Bioinformatics Core Projects/WKU-SSahi/Results/Differentially Expressed Genes/ROOT-SHOOT/")
# dir(datadir)
# setwd(datadir)
# write.table(de_TREATED_ROOT_SHOOT_20mM_ALL,file="ROOT_SHOOT_20mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_SHOOT_20mMS6Pct_ALL,file="ROOT_SHOOT_20mMS6Pct_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_SHOOT_1.25mM_ALL,file="ROOT_SHOOT_1.25mM_ALL.txt", sep="\t", col.names=T, row.names=F)
# write.table(de_TREATED_ROOT_SHOOT_P_ALL,file="ROOT_SHOOT_P_ALL.txt", sep="\t", col.names=T, row.names=F)
# 
# ###########################################
# ## NOW DO categoryCompare for Enrichment ##
# ###########################################
# 
# table_ROOT_SHOOT_20mM <- unlist(de_TREATED_ROOT_SHOOT_20mM_adj$TAIR)
# list_ROOT_SHOOT_20mM  <- list(genes=table_ROOT_SHOOT_20mM, universe=gUniverse, annotation='org.At.tair.db')
# 
# table_ROOT_SHOOT_20mMS6Pct <- unlist(de_TREATED_ROOT_SHOOT_20mMS6Pct_adj$TAIR)
# list_ROOT_SHOOT_20mMS6Pct  <- list(genes=table_ROOT_SHOOT_20mMS6Pct, universe=gUniverse, annotation='org.At.tair.db')
# 
# table_ROOT_SHOOT_P <- unlist(de_TREATED_ROOT_SHOOT_P_adj$TAIR)
# list_ROOT_SHOOT_P  <- list(genes=table_ROOT_SHOOT_P, universe=gUniverse, annotation='org.At.tair.db')
# 
# geneLists_ROOT_SHOOT <- list(T_20mM=list_ROOT_SHOOT_20mM, T_20mMS6Pct=list_ROOT_SHOOT_20mMS6Pct, T_P=list_ROOT_SHOOT_P)
# geneLists_ROOT_SHOOT <- new('ccGeneList', geneLists_ROOT_SHOOT, ccType=c('BP', 'KEGG'))
# 
# fdr(geneLists_ROOT_SHOOT) <- 0
# enrichLists_ROOT_SHOOT <- ccEnrich(geneLists_ROOT_SHOOT)
# pvalueCutoff(enrichLists_ROOT_SHOOT$BP) <- 0.001
# 
# ccOpts_ROOT_SHOOT <- new('ccOptions', listNames = names(geneLists_ROOT_SHOOT), outType='none')
# ccResults_ROOT_SHOOT <- ccCompare(enrichLists_ROOT_SHOOT, ccOpts_ROOT_SHOOT)
# cw_ROOT_SHOOT.BP <- ccOutCyt(ccResults_ROOT_SHOOT$BP, ccOpts_ROOT_SHOOT)
# breakEdges(cw_ROOT_SHOOT.BP, 0.8)
# 
# ##################################################
# ## Gene Ontology Biological Process Annotations ##
# ##################################################
# GO_BP_ROOT_SHOOT_20mM      <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT_SHOOT[[1]][[1]]), annotationDescription=goDescAnnot)
# GO_BP_ROOT_SHOOT_20mMS6Pct <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT_SHOOT[[1]][[2]]), annotationDescription=goDescAnnot)
# GO_BP_ROOT_SHOOT_P         <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT_SHOOT[[1]][[3]]), annotationDescription=goDescAnnot)
# 
# write.table(GO_BP_ROOT_SHOOT_20mM,      file="GO-BP-ROOT_SHOOT_20mM.txt", sep="\t", col.names=T, row.names=F)
# write.table(GO_BP_ROOT_SHOOT_20mMS6Pct, file="GO-BP-ROOT_SHOOT_20mMS6Pct.txt", sep="\t", col.names=T, row.names=F)
# write.table(GO_BP_ROOT_SHOOT_P,         file="GO-BP-ROOT_SHOOT_P.txt", sep="\t", col.names=T, row.names=F)
# 
# ######################
# ## KEGG ANNOTATIONS ##
# ######################
# KEGG_ROOT_SHOOT_20mM      <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT_SHOOT[[2]][[1]]), annotationDescription=KEGGDescAnnot)
# KEGG_ROOT_SHOOT_20mMS6Pct <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT_SHOOT[[2]][[2]]), annotationDescription=KEGGDescAnnot)
# KEGG_ROOT_SHOOT_P         <- getCCAnnotationTableResults(annotTable=summary(enrichLists_ROOT_SHOOT[[2]][[3]]), annotationDescription=KEGGDescAnnot)
# 
# write.table(KEGG_ROOT_SHOOT_20mM,      file="KEGG-ROOT_SHOOT_20mM.txt", sep="\t", col.names=T, row.names=F)
# write.table(KEGG_ROOT_SHOOT_20mMS6Pct, file="KEGG-ROOT_SHOOT_20mMS6Pct.txt", sep="\t", col.names=T, row.names=F)
# write.table(KEGG_ROOT_SHOOT_P,         file="KEGG-ROOT_SHOOT_P.txt", sep="\t", col.names=T, row.names=F)