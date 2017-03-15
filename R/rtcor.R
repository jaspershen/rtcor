rtcor <- function(ref.data = "RC Amide batch1 QC pos.csv",
                  mz.tolerance = 50,
                  rt.tolerance = 180,
                  max.missing = 0) {
#csv file
file <- setdiff(dir(), c(ref.data, metabolite))
file <- file[grep("csv", file)]

#read
ref.data <- read.csv(ref.data, stringsAsFactors = F)
metabolite <- read.csv(metabolite, stringsAsFactors = F)
rownames(metabolite) <- metabolite[,1]
metabolite <- metabolite[,-1]

data <- lapply(file, function(x) {read.csv(x, stringsAsFactors = F)})
names(data) <- file

ref.data1 <- ref.data[,c("name", "mzmed", "rtmed", "maxint")]
rownames(ref.data1) <- ref.data1[,1]
ref.data1 <- ref.data1[,-1]

data1 <- lapply(data, function(x) {x <- x[,c("name", "mzmed", "rtmed", "maxint")];rownames(x) <- x[,1];x <- x[,-1]})

library(MetProcesser)
ref.result <- MetProcesser::SXTMTmatch(metabolite,
                                    ref.data1,
                                    mz.tolerance = mz.tolerance,
                                    rt.tolerance = rt.tolerance)
result <-
  lapply(data1, function(x) {MetProcesser::SXTMTmatch(metabolite, x, mz.tolerance = mz.tolerance, rt.tolerance = rt.tolerance)})

ref.result1 <- data.frame(ref.result, ref.data1[ref.result[,"Index2"],"maxint"])
colnames(ref.result1)[9] <- "maxint"

index1 <- ref.result1[,"Index1"]
unique.index1 <- unique(index1)
data1.no <- setdiff(c(1:7), unique.index1)
if (length(data1.no) != 0) {
  stop(rownames(metabolite)[data1.no], "are no found in ref data!!!")
}
temp <- NULL
for (i in unique.index1) {
  temp.idx <- which(index1 %in% i)
  if (length(temp.idx) == 1) {
    temp1 <- ref.result1[temp.idx,]
    temp <- rbind(temp, temp1)} else {
      temp1 <- ref.result1[temp.idx,]
      temp1 <- temp1[which.max(temp1[,"maxint"]),]
      temp <- rbind(temp, temp1)}
}

ref.result1 <- temp

result1 <- list()
for (i in 1:length(result)) {
 result1[[i]] <- data.frame(result[[i]], data1[[i]][result[[i]][,"Index2"],"maxint"])
 colnames(result1[[i]])[9] <- "maxint"

 index2 <- result1[[i]][,"Index1"]
 unique.index2 <- unique(index2)

 data2.no <- setdiff(c(1:7), unique.index2)

 if (length(data2.no) != 0) {
   stop(rownames(metabolite)[data2.no], "are no found in ",names(result)[i],"!!!")
 }

 temp <- NULL
 for (j in unique.index2) {
   temp.idx <- which(index2 %in% j)
   if (length(temp.idx) == 1) {
     temp1 <- result1[[i]][temp.idx,]
     temp <- rbind(temp, temp1)} else {
       temp1 <- result1[[i]][temp.idx,]
       temp1 <- temp1[which.max(temp1[,"maxint"]),]
       temp <- rbind(temp, temp1)}
 }
 result1[[i]] <- temp
}

ref.rt <- ref.result1[,"rt2"]
rt <- lapply(result1, function(x) {x[,"rt2"]})

all.rt <- ref.rt
for (i in 1:length(rt)) {
  all.rt <- cbind(all.rt, rt[[i]])
}

all.rt <- cbind(metabolite[,"rt"], all.rt)
colnames(all.rt) <- c("metabolite", "reference", substr(file, 1, nchar(file) - 4))

pdf("RT distribution.pdf", width = 7, height = 7)

mycol = rainbow(n = 100)[c(1,14,23,40,59,77,95)]
par(mar = c(5,10,4,2))
plot(0,
     xlim = c(min(all.rt) - 10, max(all.rt) + 10),
     ylim = c(0, length(rt) + 3),
     xlab = "RT",
     ylab = "",
     cex.lab = 1.8,
     cex.axis=  1.5,
     pch = 19,
     col = "white",
     yaxt = "n")

axis(side = 2,
     at = c(1: (length(rt) + 2)),
     labels = colnames(all.rt),
     las = 2)
for (i in 1:nrow(all.rt)) {
  points(x = all.rt[i,],y = c(1 : (length(rt) + 2)), type = "o", pch = 19, col = mycol[i])
}

dev.off()
}


