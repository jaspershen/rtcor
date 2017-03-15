rtcor <- function(ref.data = "data1.csv",
                  cor.data = "data2.csv",
                  mz.tolerance = 20,
                  rt.tolerance = 30) {
  #csv file
  file <- dir()

  ## how to get RTQC data
  if (any(file=="RTQC.info.csv")){
    RTQC.info <- read.csv("RTQC.info.csv", stringsAsFactors = FALSE, check.names = FALSE)

    cat("RTQC you provided are:\n")
    format(RTQC.info)
    yes <- readline("Right(r) or wrong(w)?")
    if(yes=="w"){
      stop("Please change RTQC.info and try again.")
    }

    rownames(RTQC.info) <- RTQC.info[,1]
    RTQC.info <- RTQC.info[,-1]
    data1 <- read.csv(ref.data, stringsAsFactors = FALSE, check.names = FALSE)
    data2 <- read.csv(cor.data, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(data1)[which(colnames(data1) == "mzmed")] <- "mz"
    colnames(data1)[which(colnames(data1) == "rtmed") ] <- "rt"

    colnames(data2)[which(colnames(data2) == "mzmed")] <- "mz"
    colnames(data2)[which(colnames(data2) == "rtmed") ] <- "rt"

    rtqc1 <- data1[,c("name", "mz", "rt",colnames(data1)[grep("RTQC", colnames(data1))])]
    rtqc2 <- data2[,c("name", "mz", "rt",colnames(data2)[grep("RTQC", colnames(data2))])]
    ## remove peaks which are not deteced
    rtqc1 <- rtqc1[!apply(rtqc1,1,function(x) {all(x == 0)}),]
    rtqc2 <- rtqc2[!apply(rtqc2,1,function(x) {all(x == 0)}),]

    ## remove name
    rownames(rtqc1) <- as.character(rtqc1[,1])
    rownames(rtqc2) <- as.character(rtqc2[,1])
    rtqc1 <- rtqc1[,-1]
    rtqc2 <- rtqc2[,-1]

    ## find RTQC in ref data
    result1 <- SXTMTmatch(RTQC.info,
                          rtqc1,
                          mz.tolerance = mz.tolerance,
                          rt.tolerance = rt.tolerance)

    result2 <- SXTMTmatch(RTQC.info,
                          rtqc2,
                          mz.tolerance = mz.tolerance,
                          rt.tolerance = rt.tolerance)

    if (is.null(result1) | is.null(result2)) {
      stop("Aa lesat one data don't find any RTQC!")
    }

    ## matched information
    ## rtqc1 information
    index1 <- result1[,"Index1"]
    index2 <- result1[,"Index2"]
    ## one standard matched more than two peaks?
    temp <- cbind(rownames(RTQC.info)[index1],
                  result1[,c("mz error", "rt error")],
                  rtqc1[index2,-c(1,2)],
                  apply(rtqc1[index2,-c(1,2)],1,mean))
    rownames(temp) <- NULL
    colnames(temp)[1] <- "name"
    colnames(temp)[ncol(temp)] <- "int.mean"
    temp <- cbind(temp,rep(TRUE, nrow(temp)))
    colnames(temp)[ncol(temp)] <- "Which you want to remove(Change TRUR to FALSE)"

    temp <- edit(temp)
    idx <- which(temp[,ncol(temp)])

    rtqc1 <- cbind(rownames(RTQC.info)[index1[idx]],rtqc1[index2[idx],c("mz","rt")])
    colnames(rtqc1)[1] <- "name"


    ## rtqc2 information
    index1 <- result2[,"Index1"]
    index2 <- result2[,"Index2"]
    ## one standard matched more than two peaks?
    temp <- cbind(rownames(RTQC.info)[index1],
                  result2[,c("mz error", "rt error")],
                  rtqc2[index2,-c(1,2)],
                  apply(rtqc2[index2,-c(1,2)],1,mean))
    rownames(temp) <- NULL
    colnames(temp)[1] <- "name"
    colnames(temp)[ncol(temp)] <- "int.mean"
    temp <- cbind(temp,rep(TRUE, nrow(temp)))
    colnames(temp)[ncol(temp)] <- "Which you want to remove(Change TRUR to FALSE)"

    temp <- edit(temp)
    idx <- which(temp[,ncol(temp)])

    rtqc2 <- cbind(rownames(RTQC.info)[index1[idx]],rtqc2[index2[idx],c("mz","rt")])
    colnames(rtqc2)[1] <- "name"

    rtqc <- merge(rtqc1, rtqc2, by = "name")

    cat(nrow(rtqc),"out of",nrow(RTQC.info),"RTQC are detected in two data\n")
    format(rtqc)

    yes <- readline("Do you want to continute? (yes[y] or no[n])")
    if(yes=="n"){
      stop("You stop the progression.")
    }
    } else {##first station
     rtqc1 <- read.csv("RTQC1.csv", stringsAsFactors = FALSE, check.names = FALSE)
     rtqc2 <- read.csv("RTQC2.csv", stringsAsFactors = FALSE, check.names = FALSE)
     rtqc <- merge(rtqc1, rtqc2, by = "name")

     cat(nrow(rtqc),"out of",nrow(RTQC.info),"RTQC are detected in two data\n")
     format(rtqc)

     yes <- readline("Do you want to continute? (yes[y] or no[n])")
     if(yes=="n"){
       stop("You stop the progression.")
     }
    }

##correction model
  rt1 <- rtqc[,"rt.x"]
  rt2 <- rtqc[,"rt.y"]

  mse <- bestpoly(x = rt2,
                  y = rt1,
                  poly = c(1,2,3,4,5))
  idx <- which.min(mse)
  model <- lm(rt1 ~ poly(rt2,idx))

  rt1.all <- data1[,"rt"]
  rt2.all <- data2[,"rt"]

  predict.rt <- predict(object = model,
                        newdata = data.frame(rt2 = rt2.all))

  rt.cor <- rep(NA, length(predict.rt))
  index <- which(rt2.all >= rt2[1] & rt2.all <= rt2[length(rt2)])
  rt.cor[index] <- predict.rt[index]
  data2 <- data.frame(data2, rt.cor)
  write.csv(data2, "data2.new.csv")

  pdf("Experiment vs Predicted RT.pdf", width = 7, height = 7)
  par(mar = c(5,5,4,2))
  mycol = rep("black", length(rt.cor))
  mycol[index] <- "tomato"

  plot(rt2.all, predict.rt, xlim = c(0,720), ylim = c(0,720),
       xlab = "Experiment (seconds)",
       ylab = "Predicted (seconds)",
       cex.lab = 1.8,
       cex.axis = 1.5,
       col = mycol)

  abline(0, 1, lty = 2, col = "tomato")
  abline(v = rt2[1], lty = 2, col = "tomato")
  abline(v = rt2[length(rt2)], lty = 2, col = "tomato")

  plot(rt2.all, rt2.all - predict.rt,
       xlim = c(0,720),
       xlab = "Experiment (seconds)",
       ylab = "Residuals",
       cex.lab = 1.8,
       cex.axis = 1.5,
       col = mycol)
  abline(h = 0, lty = 2, col = "tomato")
  abline(v = rt2[1], lty = 2, col = "tomato")
  abline(v = rt2[length(rt2)], lty = 2, col = "tomato")
  dev.off()
}


