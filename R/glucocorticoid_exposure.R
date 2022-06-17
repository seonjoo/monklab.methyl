#' Glucocorticoid exposure score calculation for mDNA data
#'
#' Provençal, Nadine, Janine Arloth, Annamaria Cattaneo, Christoph Anacker, Nadia Cattane, Tobias Wiechmann, Simone Röh et al.
#' "Glucocorticoid exposure during hippocampal neurogenesis primes future stress response by inducing changes in DNA methylation."
#' Proceedings of the National Academy of Sciences 117, no. 38 (2020): 23280-23285.
#' @param dat : preprocessed mDNA data in the GenomicRatioSet-class format (minfi) or with the cpg site information in the rows, samples in the columns.
#'
#' @return glucocorticoid_exposure_score
#' @export
#'
#' @examples
#' data(testdata)
#' glucocorticoid_exposure(testdata)
glucocorticoid_exposure <-function(dat){
  data(Polyepigenetic.GC.exposure.weights)
  cpgindx = unlist(lapply(Polyepigenetic.GC.exposure.weights$CpG,
                          function(str){re=grep(str,row.names(dat));return(ifelse(is.null(re),NA,re))}))

  if (sum(is.na(cpgindx))>0){
    cat(sum(is.na(cpgindx)), 'cpg sites are missing. The missing cpg sites are:\n')
    print(  Polyepigenetic.GC.exposure.weights$CpG[which(is.na(cpgindx))])
    cat('The missing cpg sites are excluded from the calculation.\n')
    tmp = t(dat[cpgindx[-which(is.na(cpgindx))],]) %*%
      Polyepigenetic.GC.exposure.weights$Weights[-which(is.na(cpgindx))]
  }else{tmp = t(dat[cpgindx,]) %*% Polyepigenetic.GC.exposure.weights$Weights}
  print(tmp)

   re=data.frame(id = colnames(dat),
              glucocorticoid_exposure_score=tmp
   )

   return(re)
}
