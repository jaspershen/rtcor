rtCorrectionInHouse <- function(# ref.data = "data1.csv",
  pol = 'POS', # polarity
  inputFilePath = "data2.csv",  # result file comes from DDA Analyzer, absolute path, "/path/to/idresults.csv"
  outputFileName = NULL,
  RTCal.NEG = "RT_CAL_NEG.csv",  # RT calibration for negative
  RTCal.POS = "RT_CAL_POS.csv",  # RT calibration for positive
  method = "polyline",  # 'polyline' or 'loess'
  polyDegreeHighest = 8,  # the highest degree of polynomial regression
  degree = c(1, 2),
  translate = FALSE,
  QCFilePath = getwd(),
  outPath = getwd()) {
  #csv file1
  file <- dir()
  rtqc <- list()
  print(paste0('QC file path is: ', QCFilePath))
  print(paste0('Input path is: ', inputFilePath))
  print(paste0('Output path is: ', outPath))
  cat(paste0('This sample has ', pol, ' polarity and adjusted by ', '"', 'RT_CAL_', pol, '.csv"\n'))
  if (pol == 'POS') {
    rtqc.pos <- read.csv(file.path(QCFilePath, RTCal.POS), stringsAsFactors = FALSE, check.names = FALSE)
    rtqc <- dealQCRT(rawRTQC.info = rtqc.pos)
  }
  else if (pol == 'NEG') {
    rtqc.neg <- read.csv(file.path(QCFilePath, RTCal.NEG), stringsAsFactors = FALSE, check.names = FALSE)
    rtqc <- dealQCRT(rawRTQC.info = rtqc.neg)
  }
  else {
    stop("Please input polarity, 'POS' or 'NEG'!!!\n")
  }

  ## get RT range
  rtRangeExp <- NULL
  rtRangeStd <- NULL
  if (translate){
    rtRangeStd <- c(0.01 * 60, 11.5 * 60)
    rtRangeExp <- c(0.01 * 60, 23 * 60)
  }
  else {
    rtRangeStd <- c(0.01 * 60, 23 * 60)
    rtRangeExp <- c(0.01 * 60, 23 * 60)
  }

  ##correction model
  rtQC.std <- as.numeric(rtqc[['std']])
  rtQC.exp <- as.numeric(rtqc[['exp']])

  if (length(rtQC.std) <= polyDegreeHighest) {
    stop("'polyDegreeHighest' must be less than the points of RTQC!!!\n")
  }

  # compare and draw RT between experiment and standard
  localFileName <- file.path(outPath, paste0("Compare_RTQC_", pol, ".pdf"))
  if (!file.exists(localFileName)) {
    pdf(localFileName,
        width = 8, height = 3)
    plot(x=rtQC.std, y = rep(1, length(rtQC.std)), xlim = range(rtQC.std), ylim = c(1, 2), xlab = 'rt', ylab = "", yaxt = 'n')
    points(x = rtQC.exp, y = rep(2, length(rtQC.exp)))
    segments(x0 = rtQC.std, x1 = rtQC.exp, y0 = rep(1, length(rtQC.std)), y1 = rep(2, length(rtQC.exp)))
    grid()
    axis(side = 2, at = c(1, 2), labels = c('RT1(std)', 'RT2(exp)'))
    title(main = 'RTQC Compare')
    dev.off()
  }
  ## add the head and tail points, create this model again
  rtQC.exp <- c(rtQC.exp, rtRangeExp)
  rtQC.std <- c(rtQC.std, rtRangeStd)
  if (method == "polyline"){
    mse <- bestpoly(x = rtQC.exp,
                    y = rtQC.std,
                    pol = pol,
                    yRange = rtRangeStd,
                    xRange = rtRangeExp,
                    outPath = outPath,
                    polyDegree = polyDegreeHighest)
    idx <- which.min(mse)
    model <- lm(rtQC.std ~ poly(rtQC.exp, idx, raw = TRUE))
  }
  else {
    result <- bestloess(x = rtQC.exp,
                        y = rtQC.std,
                        span.begin = 0.5,
                        span.end = 1,
                        span.step = 0.1,
                        degree = degree,
                        path = outPath)
    idx <- which.min(result[, 3])
    model <- loess(rtQC.std ~ rtQC.exp, span = result[idx, 2], degree = result[idx, 1])
  }
  rawData <- read.csv(inputFilePath, stringsAsFactors = FALSE, check.names = FALSE)
  rtTotalExp <- as.numeric(rawData[, "rtmed"])

  predict.rt <- predict(object = model,
                        newdata = data.frame(rtQC.exp = rtTotalExp))

  rt.cor <- rep(NA, length(predict.rt))
  indMins <- which(predict.rt < min(rtRangeExp))  # lower than min RT after correction
  predict.rt[indMins] <- min(rtRangeExp)
  indMaxs <- which(predict.rt > max(rtRangeExp))  # bigger than max RT after correction
  predict.rt[indMaxs] <- max(rtRangeExp)
  # rt.cor[index] <- predict.rt[index]
  corData <- data.frame(rawData, predict.rt)
  if (is.null(outputFileName)) {
    outputFileName <- gsub('.csv', '_RTcorrected.csv', basename(inputFilePath))
  }
  write.csv(corData, file.path(outPath, outputFileName))

  localFileName <- file.path(outPath, paste0("Exp_vs_Pred_RT_", pol, ".pdf"))
  if (!file.exists(localFileName)) {
    pdf(localFileName, width = 7, height = 7)
    par(mar = c(5, 5, 4, 2))
    mycol = rep("black", length(predict.rt))
    mycol[-c(indMaxs, indMins)] <- "tomato"

    plot(rtTotalExp, predict.rt,
         xlab = "Experimental RT (s)",
         ylab = "Predicted RT (s)",
         cex.lab = 1.8,
         cex.axis = 1.5,
         col = mycol,
         pch = 19)

    abline(0, 1, lty = 2, col = "tomato")
    abline(v = min(rtQC.exp), lty = 2, col = "tomato")
    abline(v = max(rtQC.exp), lty = 2, col = "tomato")

    plot(rtTotalExp, rtTotalExp - predict.rt,
         xlab = "Experimental RT(s)",
         ylab = "Residuals (s)",
         cex.lab = 1.8,
         cex.axis = 1.5,
         col = mycol,
         pch = 19)
    abline(h = 0, lty = 2, col = "tomato")
    abline(v = rtQC.exp[1], lty = 2, col = "tomato")
    abline(v = rtQC.exp[length(rtQC.exp)], lty = 2, col = "tomato")
    dev.off()
  }
}