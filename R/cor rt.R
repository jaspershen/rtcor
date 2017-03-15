RTcorrection <- function(data1 = "data1.csv",
                         data2 = "data2.csv",
                         metabolite = "metabolite.csv",
                         mz.tolerance = 50,
                         rt.tolerance = 60) {
#读取数据
browser()
data1 <- read.csv(data1, stringsAsFactors = F)
tags1 <- data1[,c("name", "mzmed", "rtmed", "maxint")]
rownames(tags1) <- tags1[,1]
tags1 <- tags1[,-1]

data2 <- read.csv(data2, stringsAsFactors = F)
tags2 <- data2[,c("name", "mzmed", "rtmed", "maxint")]
rownames(tags2) <- tags2[,1]
tags2 <- tags2[,-1]

metabolite <- read.csv(metabolite, stringsAsFactors = F)
rownames(metabolite) <- metabolite[,1]
metabolite <- metabolite[,-1]

#在data1中寻找代谢物
result1 <- SXTMTmatch(metabolite,
                      tags1,
                      mz.tolerance = mz.tolerance,
                      rt.tolerance = rt.tolerance)

result1 <- data.frame(result1, tags1[result1[,"Index2"],"maxint"])
colnames(result1)[9] <- "maxint"

index1 <- result1[,"Index1"]
unique.index1 <- unique(index1)
tags1.no <- setdiff(c(1:7), unique.index1)
if (length(tags1.no) != 0) {
  stop(rownames(metabolite)[tags1.no], "are no found in ref data!!!")
}

temp <- NULL
for (i in unique.index1) {
  temp.idx <- which(index1 %in% i)
  if (length(temp.idx) == 1) {
    temp1 <- result1[temp.idx,]
    temp <- rbind(temp, temp1)} else {
      temp1 <- result1[temp.idx,]
      temp1 <- temp1[which.max(temp1[,"maxint"]),]
      temp <- rbind(temp, temp1)}
}

result1 <- temp

#在data2中寻找代谢物
result2 <- SXTMTmatch(metabolite,
                      tags2,
                      mz.tolerance = mz.tolerance,
                      rt.tolerance = rt.tolerance)

result2 <- data.frame(result2, tags2[result2[,"Index2"],"maxint"])
colnames(result2)[9] <- "maxint"

index1 <- result2[,"Index1"]
unique.index1 <- unique(index1)
tags2.no <- setdiff(c(1:7), unique.index1)
if (length(tags2.no) != 0) {
  stop(rownames(metabolite)[tags2.no], "are no found in ref data!!!")
}

temp <- NULL
for (i in unique.index1) {
  temp.idx <- which(index1 %in% i)
  if (length(temp.idx) == 1) {
    temp1 <- result2[temp.idx,]
    temp <- rbind(temp, temp1)} else {
      temp1 <- result2[temp.idx,]
      temp1 <- temp1[which.max(temp1[,"maxint"]),]
      temp <- rbind(temp, temp1)}
}

result2 <- temp

#开始建立模型
rt1 <- result1[,"rt2"]
rt2 <- result2[,"rt2"]

mse <- bestpoly(x = rt2,
                y = rt1,
                poly = c(1,2,3,4,5))
idx <- which.min(mse)
model <- lm(rt1 ~ poly(rt2,idx))
# library(e1071)
# model2 <- svm(rt1 ~ rt2)

rt1.all <- tags1[,"rtmed"]
rt2.all <- tags2[,"rtmed"]

predict.rt <- predict(object = model,
                      newdata = data.frame(rt2 = rt2.all))

rt.cor <- rep(NA, length(predict.rt))
index <- which(rt2.all >= rt2[1] & rt2.all <= rt2[7])
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
abline(v = rt2[7], lty = 2, col = "tomato")

plot(rt2.all, rt2.all - predict.rt,
     xlim = c(0,720),
     xlab = "Experiment (seconds)",
     ylab = "Residuals",
     cex.lab = 1.8,
     cex.axis = 1.5,
     col = mycol)
abline(h = 0, lty = 2, col = "tomato")
abline(v = rt2[1], lty = 2, col = "tomato")
abline(v = rt2[7], lty = 2, col = "tomato")
dev.off()
}


###优化参数
bestpoly <- function(x,
                     y,
                     poly = c(1, 2, 3, 4),
                     path = NULL,
                     col.list = c("palegreen", "royalblue", "firebrick1", "cyan", "tan1")

) {
  if(is.null(path)) {path <- getwd()}
  pre.all <- list()
  for (i in 1:length(poly)) {
    ##使用LOO进行参数优化
    pre <- NULL
    for(j in 1:length(x)) {
      y1 <- y[-j]
      x1 <- x[-j]
      model <- lm(y1 ~ poly(x1,poly[i]))
      pre[j] <- predict(object = model,
                        newdata = data.frame(x1 = x[j]))
    }
    pre.all[[i]] <- pre
  }

  mse <- unlist(lapply(pre.all, function(x) {sum((y -x)^2)/length(y) }))

  y.lim1 <- 0.8*min(c(y, unlist(lapply(pre.all, min))))
  y.lim2 <- 1.2*max(c(y, unlist(lapply(pre.all, max))))

  pdf(file.path(path, "MSE plot.pdf"),
      width = 7,
      height = 7)
  par(mar = c(5,5,4,2))
  plot(x, y,
       ylim = c(y.lim1, y.lim2),
       cex.lab = 1.8,
       cex.axis = 1.5,
       pch = 19,
       cex = 1.5,
       type = "o")
  abline(0, 1, lty = 2)
  for (k in 1:length(pre.all)) {
    points(x = x, y = pre.all[[k]],
           pch = 19,
           cex = 1,
           type = "o",
           col = col.list[k])
  }

  legend("topleft",
         legend = paste(paste("poly",poly),round(mse,2), sep = ":"),
         col = col.list,
         pch = 19,
         lty = 1,
         bty = "n",
         pt.cex = 1.5,
         lwd = 2,
         title = "MSE",
         cex = 1.5)
  dev.off()
  return(mse)
}


SXTMTmatch<-function(data1,data2, mz.tolerance=25, rt.tolerance = 180) {

  if (nrow(data1)==0|nrow(data2)==0) {result<-NULL;return(result)}
  mz1<-as.numeric(data1[,1])
  rt1<-as.numeric(data1[,2])

  mz2<-as.numeric(data2[,1])
  rt2<-as.numeric(data2[,2])

  result <- NULL
  cat("finished: %");cat("\n")
  for (i in 1:length(mz1)) {
    mz.error <- abs(mz1[i] - mz2)*10^6/mz1[i]
    rt.error <- abs(rt1[i] - rt2)
    j <- which(mz.error <= mz.tolerance & rt.error<=rt.tolerance)
    if (length(j) != 0) {
      result1 <- cbind(i,j,mz1[i],mz2[j],mz.error[j],rt1[i],rt2[j],rt.error[j])
      result <- rbind(result,result1)
    }

    count <- floor((length(mz1))*c(seq(0,1,0.01)))
    if (any(i==count)) {cat(ceiling (i*100/length(mz1)))
      cat(" ")}

  }
  cat("\n")
  if (is.null(result)) {cat("There are not any peak be matched\n,please change the mz or rt tolerance and try again")
    cat("\n")
  }
  else {number1 <- length(unique(result[,1]))
  number2 <- length(unique(result[,2]))
  cat(paste("There are",number1,"peaks in data1, and",number2,"peaks in data2 are matched"))
  cat("\n")
  colnames(result) <- c("Index1","Index2","mz1","mz2","mz error","rt1","rt2","rt error")
  return(result)
  }
}

