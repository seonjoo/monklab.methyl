#' Interactive QC Codes
#'
#' @param WB
#' @param idatpath
#' @param targetfilepath path of idat files
#' @param quiet F (default) do not print data import information.
#' @param sexinfo NULL (if specified, we can compare the )
#' @param minfiqc TRUE (default) run minfiqc
#' @param ewastoolqc TRUE (default) use ewastools to retrieve 17 contraol measures
#'
#' @return qc index
#' @export
#' @import minfi ewastools dplyr knitr
#' @examples
qc<-function(WB=NULL,
             idatpath=NULL,
             targetfilepath=NULL,
             targetfile=NULL,
             quiet=F,
             sexinfo=NULL,
             minfiqc=TRUE,
             ewastoolqc=TRUE,
             str='qc'

){

  cat('##############QUALITY CONTROL#################\n')


  if(minfiqc){
    cat('# minfi package\n')
    cat('# badSampleCutoff = 10.5, below the line of -log2(Meth mediation intensity) + badSampleCutoff * 2\n')

    MSet <- preprocessRaw(WB)
    out=list()
    out$qc <-  getQC(MSet)
    out$qc$predictedSex = getSex(mapToGenome(MSet))$predictedSex

    png( paste(str,'_fig1_minfi_qc.png',sep=''), width=400, height=400, res=100)
    plotQC(out$qc)
    dev.off()

  #  cat('##############SEX PREDICTION##################')
    minfi.qc = data.frame(out$qc)
    minfi.qc$truesex = sexinfo
    minfi.qc = minfi.qc %>%
      dplyr::mutate(meds = (mMed + uMed)/2) %>%
      dplyr::mutate(bad_methy = ifelse(meds < 10.5,'Fail','Pass')) %>%
      dplyr::mutate(bad_sex = ifelse(predictedSex == truesex,'Pass','Fail')) %>%
      mutate(inex=1:nrow(minfi.qc))

    if (sum(minfi.qc$bad_methy=='Fail' | minfi.qc$bad_sex=='Fail' | is.na(minfi.qc$bad_methy) | is.na(minfi.qc$bad_methy)>0, na.rm = TRUE)){
      cat(paste(sum(minfi.qc$bad_methy=='Fail' | minfi.qc$bad_sex=='Fail', na.rm = TRUE),'samples were identified as a bad sample.'))
      (minfi.qc %>% filter(bad_methy=='Fail' | bad_sex=='Fail' | is.na(bad_methy) | is.na(bad_sex)))

    }else{cat('Nothing failed.\n')}
  }


  if (ewastoolqc){
    if (quiet==F){
      ewas_meth<-invisible(read_idats( targetfilepath, quiet=quiet))
    }else{
      ewas_meth<-read_idats( targetfilepath, quiet=T)
    }
###################ILLUMINA QC Control Check#########3

# The name of the platform (450K/EPIC)
#    cat(paste('Platform:', ewas_meth$platform,'\n'))

#Looking for failed samples
    cmat = control_metrics(ewas_meth)
    targetfile$failed <- sample_failure(cmat)

    cat(paste('Number of Failed samples:',sum(targetfile$failed==TRUE),'out of ',nrow(targetfile),'samples. \n'))
 #  print(data.frame(contromat)[targetfile$failed==TRUE,])
    threshold=t(unlist(lapply(cmat, attributes))) %>% data.frame(.) %>% reshape2::melt(.)
    contromat=data.frame(id=targetfile$Basename,
                         failed = targetfile$failed,
                         failed.count = apply(as.matrix(data.frame(cmat)),1, function(x)sum(x<threshold$value, na.rm=TRUE)),
                         data.frame(cmat),
                         indx = 1:nrow(data.frame(cmat)))

    write.csv(cbind(minfi.qc,contromat), file=paste(str, '_qc_metrix.csv'), row.names=FALSE)

#' Here we see that two samples failed: 9285451058_R01C01 and 9285451058_R06C02
    if (sum(targetfile$failed==TRUE)>0){

  # below threshold will be considered as problematic.

 #paste(names(contromat)[-1], contromat_count)

      failedindx = which(targetfile$failed==TRUE)
      cat(paste('Failed Sample Index: ',failedindx, '\n'))

      threshold$variable=gsub('.threshold','',threshold$variable)

      cbind(t(contromat %>% filter(failed==TRUE)),
                      c('threshold', '',' ',threshold$value))

    }
  }
}
