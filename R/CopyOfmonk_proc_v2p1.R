#' DMN methylation preprocessing
#' Epi 450K Data Pre-processing is first established on June 10, 2021 By Pamela Scorza
#' Version 2 monk preprocessing Last update: 2026 09 02: In the process of revision,
#' preprocessing pipeline updated is required. We updated as following:
#'
#' 1. adding a bead count threshold, and for both detP and bead count,
#' removing probes that fail either in >5% of samples. (in previous version, we excluded )
#'
#' 2. add Genomic-based probe filtering

#' @param WB methylation file imported by minfi package
#' @param idatpath  idat file path
#' @param targetfile sample information file. This must include Basename field that include idat_id
#' @param probthresh 0.01 (default) p-value threshold to filter out significantly methylation CPG sites.
#' @param bead_min_count 3 (default)
#' @param fail_fraction 0.05 (default)
#' @param snp_maf_dropLoci 0 (default)
#' @param crosshyb_maf 0.05 (default)
#' @param rmXY FALSE (default) remove XY chromosomes
#' @param sample_fail_fraction 0.05 (default) In probe-level filtering, probes where qc failire happenes greater than sample_fail_fraction,
#'will be removed
#' @param median_intensity_cutoff  10.5 (default) for removing sample.
#' @param outfilename monkproc_v2 (devault) output file name
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
#' #WB <- read.metharray.exp(base=idatPath, targets=targetfile, verbose=T, extended=TRUE) # read the idat file one by one

#' #monk_proc_v2(idatpath,targetfile,probthresh=0.01,outfilename='monkproc_v2_test')


monk_proc_v2 <- function(WB,
                         idatpath = NULL,
                         targetfile = NULL,
                         probthresh = 0.01,
                         bead_min_count = 3,
                         fail_fraction = 0.05,
                         snp_maf_dropLoci = 0,
                         crosshyb_maf = 0.05,
                         rmXY = FALSE,
                         sample_fail_fraction = 0.05,
                         median_intensity_cutoff = 10.5,
                         outfilename = 'monkproc_v2'
){

  ######################################################################
  ## Detect array type
  ######################################################################
  array_type <- minfi::annotation(WB)["array"]
  is_epic <- grepl("EPIC", array_type, ignore.case = TRUE)

  cat("Start Mproc Preprocessing v2 (2026 09 02 revision).\n")
  cat(paste('Detected array type:', array_type, '\n'))
  cat(paste('Number of imported methylation samples:', ncol(WB), '\n'))

  ######################################################################
  ## Sample-level QC — computed entirely from WB, no re-reading of IDATs
  ######################################################################
  cat('## Sample-level QC\n')

  meta <- targetfile

  # Assumes meta row order matches colnames(WB) — verify against your target file.

  ## --- 17 BeadArray control metrics, extracted from WB's control probes ---
  ewas_meth<-read_idats( paste(idatpath, '/',targetfile$Basename,sep=''), quiet=T)

  cmat = control_metrics(ewas_meth)
  threshold=t(unlist(lapply(cmat, attributes))) %>% data.frame(.) %>% reshape2::melt(.)

  contromat=data.frame(id=targetfile$Basename,
#                       failed = targetfile$failed,
                       failed.count = apply(as.matrix(data.frame(cmat)),1, function(x)sum(x<threshold$value, na.rm=TRUE)),
                       data.frame(cmat),
                       indx = 1:nrow(data.frame(cmat)))

  ## --- Median methylation intensity ---
  qc <- getQC(preprocessRaw(WB))

  ## --- Predicted sex ---
  sex <- getSex(mapToGenome(preprocessRaw(WB)))

  ## --- Combine into one sample QC table and export ---
  sampleQC <- data.frame(
    contromat,
    mMed           = qc$mMed,
    uMed           = qc$uMed,
    medianInt_fail = qc$mMed < median_intensity_cutoff | qc$uMed < median_intensity_cutoff,
    predictedSex   = sex$predictedSex,
    sex_pred_fail = targetfile$baby_sex != sex$predictedSex
  )
  sampleQC$overall_fail <- with(sampleQC, failed.count>2 | sex_pred_fail | medianInt_fail)

  write.csv(sampleQC,
            paste0(outfilename, "_sampleQC_", ifelse(is_epic, "EPIC", "450k"), ".csv"),
            row.names = FALSE)

  WB <- WB[, !sampleQC$overall_fail]
  cat(paste(sum(sampleQC$overall_fail), 'sample(s) excluded;', ncol(WB), 'remain\n'))

  rm(ctrlInfo, grn, red, detP, qc, sex); gc()

  ######################################################################
  ## Normalization
  ######################################################################
  WB.noob <- preprocessNoob(WB, dyeMethod = "single")

  cat('# Distribution of beta-values: before and after noob normalization\n')
  densityPlot(WB, main = paste("density plots before/after preprocessing —", array_type),
              pal = "#440154FF", ylim = c(0, 4.5))
  densityPlot(WB.noob, add = FALSE, pal = "#FDE725FF")
  legend("topleft", c("Noob", "Raw"), lty = c(1, 1), title = "Normalization",
         bty = 'n', cex = 1.3, col = c("#FDE725FF", "#440154FF"))

  ######################################################################
  ## Probe-level QC: bead count + detection p-value
  ######################################################################
  cat('## Remove bad/failed probes (bead count, detection p-value)\n')

  filter_good_beadcount_probes <- function(rgset, probes = featureNames(rgset),
                                           bead_min_count = 3, bead_fail_fraction = 0.05) {
    nbeads <- wateRmelon::beadcount(rgset)
    probes <- intersect(probes, rownames(nbeads))
    nbeads <- nbeads[probes, , drop = FALSE]
    keep <- rowMeans(is.na(nbeads) | nbeads < bead_min_count, na.rm = TRUE) <= bead_fail_fraction
    names(keep)[keep]
  }

  filter_good_detection_probes <- function(rgset, probes,
                                           detection_p_threshold = 0.01, detection_p_fail_fraction = 0.05) {
    det_p <- minfi::detectionP(rgset)
    probes <- intersect(probes, rownames(det_p))
    det_p <- det_p[probes, , drop = FALSE]
    keep <- rowMeans(is.na(det_p) | det_p > detection_p_threshold, na.rm = TRUE) <= detection_p_fail_fraction
    names(keep)[keep]
  }

  all_probes <- rownames(getAnnotation(WB))
  good_beads <- filter_good_beadcount_probes(WB, all_probes, bead_min_count, fail_fraction)
  good_detection <- filter_good_detection_probes(WB, all_probes, probthresh, fail_fraction)

  before <- nrow(WB.noob)
  WB.noob <- WB.noob[rownames(WB.noob) %in% good_beads, ]
  WB.noob <- WB.noob[rownames(WB.noob) %in% good_detection, ]

  cat(paste('[', array_type, '] Started with', before, 'probes;', nrow(WB.noob),
            'remain after bead count/detection p-value filtering —',
            before - nrow(WB.noob), 'removed (',
            round(100 - nrow(WB.noob) / before * 100, 2), '%)\n'))

  ######################################################################
  ## Genomic mapping + SNP removal
  ######################################################################
  cat('## Map to genome and remove SNP-affected probes\n')

  WB.noob <- mapToGenome(WB.noob)
  WB.noob <- minfi::dropLociWithSnps(WB.noob, snps = c("CpG", "SBE"), maf = snp_maf_dropLoci, snpAnno = NULL)

  ######################################################################
  ## Cross-reactive probe removal — array-specific list
  ######################################################################
  cat('## Remove cross-reactive/polymorphic probes\n')

  before <- nrow(WB.noob)
  betas_tmp <- getBeta(WB.noob)

  if (!is_epic) {
    cat('Using wateRmelon::rmSNPandCH (450k cross-reactive list, Chen et al. 2013)\n')
    betas_tmp <- rmSNPandCH(betas_tmp, mafcut = crosshyb_maf, and = TRUE,
                            rmcrosshyb = TRUE, rmXY = rmXY)
  } else {
    cat('Using sesame EPIC masking (Zhou et al. 2017 / Pidsley et al. 2016 lineage)\n')
    if (!requireNamespace("sesameData", quietly = TRUE) || !requireNamespace("sesame", quietly = TRUE)) {
      stop("EPIC cross-reactive filtering requires the 'sesame' and 'sesameData' packages.")
    }
    mask_probes <- sesameData::sesameDataGet("EPIC.probeInfo")$mask
    keep_probes <- setdiff(rownames(betas_tmp), mask_probes)
    betas_tmp <- betas_tmp[keep_probes, , drop = FALSE]

    if (rmXY) {
      ann_epic <- getAnnotation(WB.noob)
      xy_probes <- rownames(ann_epic)[ann_epic$chr %in% c("chrX", "chrY")]
      betas_tmp <- betas_tmp[!(rownames(betas_tmp) %in% xy_probes), , drop = FALSE]
    }
  }
  # rmXY = FALSE retains sex-chromosome probes by default — confirm this is
  # the intended choice for your analysis.

  WB.noob <- WB.noob[rownames(WB.noob) %in% rownames(betas_tmp), ]
  cat(paste('[', array_type, ']', before - nrow(WB.noob),
            'cross-reactive/polymorphic probes removed;', nrow(WB.noob), 'probes remain\n'))
  rm(betas_tmp); gc()

  ######################################################################
  ## Probe-type bias correction (RCP)
  ######################################################################
  cat('# Probe-type bias correction: RCP (Niu et al., Bioinformatics 2016)\n')
  MSet_tmp <- MethylSet(Meth = getMeth(WB.noob), Unmeth = getUnmeth(WB.noob),
                        colData = colData(WB.noob), annotation = annotation(WB.noob))
  betas.rcp <- rcp(MSet_tmp)

  ######################################################################
  ## Batch effects
  ######################################################################
  cat('## Batch effects\n')

  pData(WB.noob)$Sentrix_ID <- sapply(colnames(WB.noob), function(ss) strsplit(ss, '_')[[1]][1])
  pData(WB.noob)$array_row  <- substring(sapply(colnames(WB.noob), function(ss) strsplit(ss, '_')[[1]][2]), 1, 3)
  pData(WB.noob)$array_col  <- substring(sapply(colnames(WB.noob), function(ss) strsplit(ss, '_')[[1]][2]), 4, 7)
  pData(WB.noob)$array_rowcol <- substring(sapply(colnames(WB.noob), function(ss) strsplit(ss, '_')[[1]][2]), 1, 7)

  PCobject <- prcomp(t(betas.rcp), retx = TRUE, center = TRUE, scale. = TRUE)
  PCs <- PCobject$x

  cat('# Association between PC1 and array position / chip\n')
  # Due to the way the sample was run, we don't have sufficient numbersof chips..
  #try(print(oneway.test(PCs[, 1] ~ as.factor(pData(WB.noob)$Sentrix_ID))))

  try(print(oneway.test(PCs[, 1] ~ as.factor(pData(WB.noob)$array_row))))
  try(print(oneway.test(PCs[, 1] ~ as.factor(pData(WB.noob)$array_col))))
  print(boxplot(PCs[, 1] ~ pData(WB.noob)$array_row, ylab = "PC1", las = 2,
                main = paste("Row —", array_type), col = rainbow(8)))

  # NOTE: before ComBat, confirm array position/chip is not confounded with
  # the exposure/outcome variable(s) of interest for this analysis.

  Mvals <- log2(betas.rcp) - log2(1 - betas.rcp)
  Mvals.ComBat <- ComBat(Mvals, batch = pData(WB.noob)$array_rowcol)
#  Mvals.ComBat <- ComBat(Mvals.ComBat, batch = pData(WB.noob)$Sentrix_ID)
  betas.clean <- 2^Mvals.ComBat / (1 + 2^Mvals.ComBat)

  PC_post <- prcomp(t(betas.clean), retx = TRUE, center = TRUE, scale. = TRUE)
  par(mfrow = c(1, 2))
  boxplot(PC_post$x[, 1] ~ pData(WB.noob)$array_row, ylab = "PC1", las = 2,
          main = "Row (post-correction)", col = rainbow(8))
  boxplot(PC_post$x[, 1] ~ pData(WB.noob)$array_col, ylab = "PC1", las = 2,
          main = "Column (post-correction)", col = rainbow(8))

# Due to the way the sample was run, we don't have sufficient numbersof chips..
#  par(mfrow = c(1, 1))
#  print(boxplot(PC_post$x[, 1] ~ pData(WB.noob)$Sentrix_ID, ylab = "PC1", las = 2,
#                main = "Chip (post-correction)", col = rainbow(8)))
#  try(print(oneway.test(PC_post$x[, 1] ~ as.factor(pData(WB.noob)$Sentrix_ID))))

  rm(PCs, PCobject, PC_post, Mvals, Mvals.ComBat); gc()

  ######################################################################
  ## Save
  ######################################################################
  outfilename_full <- paste0(outfilename, "_", ifelse(is_epic, "EPIC", "450k"), ".RData")
  save(betas.clean, targetfile, array_type, file = outfilename_full)
  cat(paste('Cleaned file saved as', outfilename_full, '\n'))
  cat('Done!\n')

  invisible(betas.clean)
}
