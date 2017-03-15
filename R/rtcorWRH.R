rtcorWRH <- function(ref.data = "data1.csv",
                     cor.data = "data2.csv") {
  #csv file
  file <- dir()

    RTQC.info1 <- read.csv("RTQC.info1.csv", stringsAsFactors = FALSE, check.names = FALSE)
    RTQC.info2 <- read.csv("RTQC.info2.csv", stringsAsFactors = FALSE, check.names = FALSE)

    rtqc.name1 <- unique(RTQC.info1[,1])
    rtqc.name2 <- unique(RTQC.info2[,1])

    rtqc1.ID <- NULL
    rtqc1.rt <- NULL
    for(i in 1:length(rtqc.name1)) {
      temp.idx <- which(RTQC.info1[,1] == rtqc.name1[i])
      rtqc1.ID[i] <- rtqc.name1[i]
      rtqc1.rt[i] <- mean(RTQC.info1[temp.idx,2])
    }

    rtqc1 <- cbind(rtqc1.ID, rtqc1.rt)
    rtqc1 <- rtqc1[order(rtqc1[,1]),]

    rtqc2.ID <- NULL
    rtqc2.rt <- NULL
    for(i in 1:length(rtqc.name2)) {
      temp.idx <- which(RTQC.info2[,1] == rtqc.name2[i])
      rtqc2.ID[i] <- rtqc.name2[i]
      rtqc2.rt[i] <- mean(RTQC.info2[temp.idx,2])
    }

    rtqc2 <- cbind(rtqc2.ID, rtqc2.rt)
    rtqc2 <- rtqc2[order(rtqc2[,1]),]


  ##correction model
  rt1 <- as.numeric(rtqc1[,2])
  rt2 <- as.numeric(rtqc2[,2])

  mse <- bestpoly(x = rt2,
                  y = rt1,
                  poly = c(1,2,3,4,5))
  idx <- which.min(mse)
  model <- lm(rt1 ~ poly(rt2,idx))

  data1 <- read.csv(ref.data, stringsAsFactors = FALSE, check.names = FALSE)
  data2 <- read.csv(cor.data, stringsAsFactors = FALSE, check.names = FALSE)

  rt1.all <- as.numeric(data1[,"rt"])
  rt2.all <- as.numeric(data2[,"rt"])

  predict.rt <- predict(object = model,
                        newdata = data.frame(rt2 = rt2.all))

  rt.cor <- rep(NA, length(predict.rt))
  index <- which(rt2.all >= min(rt2) & rt2.all <= max(rt2))
  rt.cor[index] <- predict.rt[index]
  data2 <- data.frame(data2, rt.cor)
  write.csv(data2, "data2.new.csv")

  pdf("Experiment vs Predicted RT.pdf", width = 7, height = 7)
  par(mar = c(5,5,4,2))
  mycol = rep("black", length(rt.cor))
  mycol[index] <- "tomato"

  plot(rt2.all, predict.rt,
       xlab = "Experiment (seconds)",
       ylab = "Predicted (seconds)",
       cex.lab = 1.8,
       cex.axis = 1.5,
       col = mycol,
       pch = 19)

  abline(0, 1, lty = 2, col = "tomato")
  abline(v = min(rt2), lty = 2, col = "tomato")
  abline(v = max(rt2), lty = 2, col = "tomato")

  plot(rt2.all, rt2.all - predict.rt,
       xlab = "Experiment (seconds)",
       ylab = "Residuals",
       cex.lab = 1.8,
       cex.axis = 1.5,
       col = mycol,
       pch = 19)
  abline(h = 0, lty = 2, col = "tomato")
  abline(v = rt2[1], lty = 2, col = "tomato")
  abline(v = rt2[length(rt2)], lty = 2, col = "tomato")
  dev.off()
}


