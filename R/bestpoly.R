###parameters optimization
bestpoly <- function(x,
                     y,
                     poly = c(1, 2, 3, 4),
                     path = NULL,
                     col.list = c("palegreen", "royalblue", "firebrick1", "cyan", "tan1")

) {
  if(is.null(path)) {path <- getwd()}
  pre.all <- list()
  for (i in 1:length(poly)) {
    ##LOO
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