#' DMN methylation preprocessing
#' Epi 450K Data Pre-processing is first established on June 10, 2021 By Pamela Scorza
#' Version 1 monk preprocessing Last update: 2022 04 08
#' @param WB methylation file imported by minfi package
#' @param idatpath idat file path
#' @param targetfile sample information file. This must include Basename field that include idat_id
#' @param probthresh 0.01 (default) p-value threshold to filter out significantly methylation CPG sites.
#' @param outfilename monkproc_v1 (devault) output file name
#'
#' @return The resulting beta
#' @export
#' @import minfi ENmix sva DMRcate dplyr
#' @examples
#' #library(monklab.methyl)
#' #library(dplyr)
#' #idatPath <- "/Volumes/ALSPAC/EPI/EPI_methylation_rawdata_froam_ben_2019_12_20/idat" # path of the folder
#' #targets <- readxl::read_excel("/Volumes/ALSPAC/EPI/EPI_methylation_rawdata_froam_ben_2019_12_20/monk_sample_450K_placenta.xlsx")
#' #targetfile=targets[1:5,c(1,2,3,4,11,12,13,14,15,17)]
#' #targetfile$Basename <- targetfile$idat_id#paste0(targets$sentrixbarcode, "_", targets$samplesection) # name of the files
#' #WB <- read.metharray.exp(base=idatPath, targets=targetfile, verbose=T) # read the idat file one by one

#' #monk_proc_v1(idatpath,targetfile,probthresh=0.01,outfilename='monkproc_v1_test')
monk_proc_v1<-function(WB,
                       idatpath=NULL,
                       targetfile=NULL,
                       probthresh=0.01, # p.detect threshold of the bad probes
                       outfilename='monkproc_v1'
                    ){


  ######################################################################
  ###################Reading in idat files##############################
  ######################################################################

  cat(paste('Number of the imported methylation samples: ',ncol(WB), '\n'))

  WB.noob <- preprocessNoob(WB)

  cat('# Distribution of beta-values: before and after noob normalization\n')

  #############FIGURE 1######################
  #' Distribution of beta-values: before and after noob normalization
  #+ fig.width=8, fig.height=6, dpi=300
  #+
  densityPlot(WB, main = "density plots before and after preprocessing", pal="#440154FF", ylim=c(0,4.5))
  densityPlot(WB.noob, add = F, pal = "#FDE725FF")
  # Add legend
  legend("topleft", c("Noob","Raw"),
         lty=c(1,1), title="Normalization",
         bty='n', cex=1.3, col=c("#FDE725FF","#440154FF"))
  #' notice the blue density traces (raw) are more spread out; background correction brings them together
  #'

  cat('## Remove Bad/Failed Probe')

  # We want to drop probes with intensity that is not significantly above background signal (from negative control probes)

  # probe failures due to low intensities

  # So, high p-values can be interpreted as "not different from background" (i.e. "no bueno")


  detect.p <- minfi::detectionP(WB, type = "m+u")
  #detect.p2 <- ewastools::detectionP(ewas_meth) #this is the other way to do detection p-value, using ewastools. Really, it's whatever you feel like doing. Note that you need the ewastools version of the file, not the WB file
  #save.image('/Volumes/ALSPAC/EPI/Monk_DNA_Methylation_data_Sameera_2022_01_20/epi_epic_placenta_preprocessed_2022_02_24.Rdata')

  cat(paste(max(colSums(detect.p > 0.01))),'probes undetected per sample if we use 10^-2 cutoff\n')

  cat(paste(max(colSums(detect.p > 0.000001))),'probes undetected per sample if we use 10^-6 cutoff\n')

  # Barplot of the mean detection p-value by sample
  #############FIGURE 2######################
  #'Barplot of the mean detection p-value by sample
  #'This tells us something (??) about the sample
  barplot(colMeans(detect.p), col=rainbow(dim(detect.p)[2]), las=2, cex.names=0.7, main="Mean detection P by sample",cex.axis=0.8, ylim=c(0,7e-4))
  #Saved this as 40*20

  # Restrict data to good probes only:

  detect.p[detect.p > probthresh] <- NA
  detect.p <- na.omit(detect.p)
  intersect <- intersect(rownames(getAnnotation(WB)), rownames(detect.p))


  cat(paste('Started with', nrow(getAnnotation(WB)), 'probes and', length(intersect),
            ' remaining if we use', probthresh,'as cutoff. Removed',nrow(getAnnotation(WB))-length(intersect),'probes',
            round(100 - length(intersect)/nrow(getAnnotation(WB))*100,2), '%\n'))


  # Filter bad probes from our methylset

  #nrow(WB.noob)
  WB.noob <- WB.noob[rownames(getAnnotation(WB.noob)) %in% intersect,]
  nrow(WB.noob)
  rm(intersect, detect.p); gc() # cleanup

  cat('######################################################################\n')
  cat('#####################Probe type adjustment############################\n')
  cat('######################################################################\n')

  cat('# Need to adjust for probe-type bias Infinium I (type I) and Infinium II (type II) probes\n')
  cat('# RCP with EnMix: Regression on Correlated Probes [Niu et al. Bioinformatics 2016](http://www.ncbi.nlm.nih.gov/pubmed/27153672)\n')
  betas.rcp <- rcp(WB.noob)
  cat('# note that this package takes beta values out of the minfi object - result is a matrix\n')


  ## Annotation of Infinium type for each probe (I vs II)
  typeI <-   minfi::getProbeInfo(WB.noob,type="I")$Name
  typeII <-  minfi::getProbeInfo(WB.noob,type="II")$Name
  onetwo <- rep(1, nrow(betas.rcp))
  onetwo[rownames(betas.rcp) %in% typeII] <- 2


  cat('# Probe-type bias adjustment before and after RCP\n')
  #+ fig.width=15, fig.height=7, dpi=300
  par(mfrow=c(1,2)) # Side-by-side density distributions
  densityPlot(WB.noob[rownames(getAnnotation(WB.noob)) %in% typeI,],pal = "#FDE725FF",main='Beta density',ylim=c(0,6.5))
  densityPlot(WB.noob[rownames(getAnnotation(WB.noob)) %in% typeII,],add = F, pal = "#440154FF")
  densityPlot(betas.rcp[rownames(getAnnotation(WB.noob)) %in% typeI,],pal = "#FDE725FF",main='Beta density probe-type adjusted',ylim=c(0,6.5))
  densityPlot(betas.rcp[rownames(getAnnotation(WB.noob)) %in% typeII,],add = F, pal = "#440154FF")
  legend("center", c("Infinium I","Infinium II"),
         lty=c(1,1), title="Infinium type",
         bty='n',col=c("#FDE725FF","#440154FF"))
  #' notice that the type I and II peaks are more closely aligned after rcp adjustment
  #' (particularly in the higher peak)
  rm(onetwo, typeI, typeII)

  cat('## Batch effects\n')

  ######################################################################
  ###########################Batch effects##############################
  ######################################################################

  #' ## Batch effects
  #' This can create batch effects (technical variation) with different intensities by position (row effect).
  #' Other commonly observed batch effects include bisulfite processing plate, chip, and processing date.

  #' We will test 3 potential factors
  pData(WB.noob)$Sentrix_ID <- unlist(lapply(colnames(WB.noob), function(ss)strsplit(ss,'_')[[1]][1]))
  pData(WB.noob)$array_row <- substring(unlist(lapply(colnames(WB.noob), function(ss)strsplit(ss,'_')[[1]][2])),1,3) #Chip row
  pData(WB.noob)$array_col <- substring(unlist(lapply(colnames(WB.noob), function(ss)strsplit(ss,'_')[[1]][2])),4,7) #Chip row #Chip position
  pData(WB.noob)$array_rowcol <- substring(unlist(lapply(colnames(WB.noob), function(ss)strsplit(ss,'_')[[1]][2])),1,7)

    #and just chip (Sentrix_ID)

  #' ## Principal Component Analysis (PCA)
  #' Calculate major sources of variability of DNA methylation using PCA
  PCobject <- prcomp(t(betas.rcp), retx = T, center = T, scale. = T)
  #' Extract the Principal Components from SVD
  PCs <- PCobject$x

  cat('# Is the major source of variability associated with position on chip?')
#  summary(lm(PCs[, 1] ~ pData(WB.noob)$array_row)) #some effects w/row
#  summary(lm(PCs[, 1] ~ pData(WB.noob)$array_col)) #no effects w/col
  print(oneway.test(PCs[, 1] ~ as.factor(pData(WB.noob)$Sentrix_ID))) #
  print(oneway.test(PCs[, 1] ~ as.factor(pData(WB.noob)$array_row))) #
  try(print(oneway.test(PCs[, 1] ~ as.factor(pData(WB.noob)$array_col)))) #

  print(boxplot(PCs[, 1] ~ pData(WB.noob)$array_row, ylab = "PC1",las=2, main="Row",col=rainbow(8)))

  #' ComBat eBayes adjustment using a known variable of interest (here we use row)
  Mvals <- log2(betas.rcp)-log2(1-betas.rcp)

  combatoption=menu(c("row", "col", 'both', 'none'), title="Which batch effect correction?")
  if (combatoption==1){
    cat('only correct for rows.\n')
    Mvals.ComBat <- ComBat(Mvals, batch = pData(WB.noob)$array_row)

  }
  if (combatoption==2){
    cat('only correct for columns.\n')
    Mvals.ComBat <- ComBat(Mvals, batch = pData(WB.noob)$array_col)

  }
  if (combatoption==3){
    cat('both columns and raws\n')
    Mvals.ComBat <- ComBat(Mvals, batch = pData(WB.noob)$array_rowcol)

  }
  if (combatoption==4){
    cat('No position effect correction.\n')
    Mvals.ComBat <- Mvals
  }

  betas.rcp <- 2^Mvals.ComBat/(1+2^Mvals.ComBat) # Convert M-values back to beta-values

  #' PCA after removing batch effects
  PC_post_correction <- prcomp(t(betas.rcp), retx = T, center = T, scale. = T)
  PC1_post_correction <- PC_post_correction$x

  par(mfrow = c(1, 2))
  boxplot(PC1_post_correction[, 1] ~ pData(WB.noob)$array_row, ylab = "PC1",las=2, main="Row",col=rainbow(8))
  boxplot(PC1_post_correction[, 1] ~ pData(WB.noob)$array_col, ylab = "PC1",las=2, main="Column",col=rainbow(8))

  #'By Beadchip
  oneway.test(PC1_post_correction[, 1] ~ as.factor(pData(WB.noob)$Sentrix_ID)) #

  table(pData(WB.noob)$Sentrix_ID)

  par(mfrow = c(1, 1))
  boxplot(PC1_post_correction[, 1] ~ pData(WB.noob)$Sentrix_ID, ylab = "PC1",las=2, main="Chip",col=rainbow(8))

  chipmenu=menu(c("no","yes"), title="Batch correction for chips?")
  if (chipmenu==1){
    cat('No chip correction.\n')
  }
  if (chipmenu==2){
    cat('Chip correction.\n')
    Mvals.ComBat <- ComBat(Mvals.ComBat, batch = pData(WB.noob)$Sentrix_ID)
    betas.rcp <- 2^Mvals.ComBat/(1+2^Mvals.ComBat) # Convert M-values back to beta-values
  }

  rm(PCs, PCobject, PC_post_correction, PC1_post_correction); gc()


  cat('######################################################################\n')
  cat('########################Cross-Hybridizing#############################\n')
  cat('######################################################################\n')

  betas.clean <- rmSNPandCH(betas.rcp,  mafcut = 0.05, and = TRUE, rmcrosshyb = TRUE, rmXY= FALSE)
  Mvals.clean <- log2(betas.clean)-log2(1-betas.clean)

  #saveRDS(betas.clean, file = "tcell_epic_n152_veronica_pamela_clean_betas.rds")
  #saveRDS(Mvals.clean, file = "tcell_epic_n152_veronica_pamela_Clean_mvals.rds")


  save(betas.clean, targetfile, file = outfilename)


  cat(paste('The cleaned file was saved as', outfilename),'\n')


  cat('Done!')
}
