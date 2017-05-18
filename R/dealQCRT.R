dealQCRT <- function(rawRTQC.info) {
  rtqc <- list()
  rtqc[['name']] <- unique(rawRTQC.info$LABID)
  # convert mintues to seconds
  if(all(rawRTQC.info$RT_EXP < 240)) {
    rawRTQC.info$RT_EXP <- rawRTQC.info$RT_EXP * 60
  }
  if(all(rawRTQC.info$RT_STD < 240)) {
    rawRTQC.info$RT_STD <- rawRTQC.info$RT_STD * 60
  }
  
  for(i in 1:length(rtqc[['name']])) {
    temp.idx <- which(rawRTQC.info$LABID == rtqc[['name']][i])
    rtqc[['labID']][i] <- rtqc[['name']][i]
    # may more than one RT
    rtqc[['exp']][i] <- mean(rawRTQC.info$RT_EXP[temp.idx]) # experimental RT
    rtqc[['std']][i] <- mean(rawRTQC.info$RT_STD[temp.idx]) # standard RT
  }
  return(rtqc)
}
