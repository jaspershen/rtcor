###parameters optimization
bestpoly <- function(x,
                     y,
                     pol = '',
                     polyDegree = 4,
                     outPath = NULL,
                     yRange = c(0.05 * 60, 23 * 60),
                     xRange = c(0.05 * 60, 23 * 60),
                     col.list = c("palegreen", "royalblue", "firebrick1", "cyan", "tan1")

) {
  if(is.null(outPath)) {outPath <- getwd()}
  pre.all <- list()
  # browser()
  for (i in seq(polyDegree)) {
    ##Leave One Out
    pre <- NULL
    for (j in 1:(length(x) - 2)) {
      y1 <- y[-j]
      x1 <- x[-j]
      model <- lm(y1 ~ poly(x1, seq(polyDegree)[i]))
      pre[j] <- predict(object = model,
                        newdata = data.frame(x1 = x[j]))
    }
    pre.all[[i]] <- pre
  }

  mse <- unlist(lapply(pre.all,
                       function(x) {sum((head(y, -2) -x)^2)/length(head(y, -2))}))
  # y axis scope
  x.lim1 <- 0.8*0.05*60
  x.lim2 <- 1.2*23*60
  y.lim1 <- 0.8*min(c(y, unlist(lapply(pre.all, min))))
  y.lim2 <- 1.2*max(c(y, unlist(lapply(pre.all, max))))

  ## Predict Model MSE Plot
  localFileName <- file.path(outPath, paste0("Pred_MSE_plot_", pol, ".pdf"))
  if (!file.exists(localFileName)) {
    pdf(localFileName,
        width = 7,
        height = 7)
    par(mar = c(5,5,4,2))
    plot(x, y,
         ylim = c(y.lim1, y.lim2),
         xlim = c(x.lim1, x.lim2),
         cex.lab = 1.8,
         cex.axis = 1.5,
         pch = 19,
         cex = 1.5,
         type = "o")
    abline(0, 1, lty = 2)
    xOutLastTow <- head(x, -2)
    for (k in 1:length(pre.all)) {
      # browser()
      points(x = xOutLastTow, y = pre.all[[k]],
             pch = 19,
             cex = 1,
             type = "o",
             col = col.list[k])
    }
    legend("topleft",
           legend = paste(paste("poly", seq(polyDegree)), round(mse, 2), sep = ": "),
           col = col.list,
           pch = 19,
           lty = 1,
           bty = "n",
           pt.cex = 1.5,
           lwd = 2,
           title = "Leave One Out MSE",
           cex = 1.5)
    dev.off()
  }
  ## selected model MSE
  ## add the head and tail points
  localFileName <- file.path(outPath, paste0("Selected_Model_MSE_plot_", pol, ".pdf"))
  if (!file.exists(localFileName)) {
    pdf(localFileName,
        width = 7,
        height = 7)
    par(mar = c(5,5,4,2))
    plot(x, y,
         ylim = c(y.lim1, y.lim2),
         xlim = c(x.lim1, x.lim2),
         cex.lab = 1.8,
         cex.axis = 1.5,
         pch = 19,
         cex = 1.2,
         type = "o")
    abline(0, 1, lty = 2)
    idx <- which.min(mse)
    bestModel <- lm(y ~ poly(x, idx, raw = TRUE))
    lines(seq(x.lim1, x.lim2, by = 0.1),
          predict(bestModel, data.frame(x = seq(x.lim1, x.lim2, by = 0.1))),
          col=col.list[idx],
          cex = 1.5)
    labList <- list()
    for (i in seq(idx + 1)) {
      labList[i] = paste0(sprintf("%0.4g", bestModel$coefficients[[i]]),
                          paste0('x^', i - 1))
    }
    browser()
    labelString <- paste0('y = ', paste(rev(labList), collapse = ' + '))
    title(main = "Selected Model")
    text(200, 10, adj = c(0, 0), labels = labelString)
    legend("topleft",
           legend = paste(paste("poly", idx),
                          round(unlist(mse[idx]), 2), sep = ": "),
           col = col.list,
           pch = 19,
           lty = 1,
           bty = "n",
           pt.cex = 1.5,
           lwd = 2,
           title = 'MSE',
           cex = 1.5)
    dev.off()
  }
  return(mse)
}