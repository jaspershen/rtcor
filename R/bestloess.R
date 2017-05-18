###parameters optimization
bestloess <- function(x,
                      y,
                      span.begin = 0.5,
                      span.end = 1,
                      span.step = 0.1,
                      degree = c(1,2),
                      col.list = rainbow(20),
                      path = NULL) {
  if(is.null(path)) {path <- getwd()}
  span <- seq(span.begin, span.end, span.step)
  para <- NULL
  for(i in 1:length(degree)){
    para <- rbind(para, cbind(degree[i], span))
  }
  colnames(para) <- c("degree", "span")
  
  y <- y[order(x)]
  x <- sort(x)
  
  pre.all <- list()
  
  for(i in 1:nrow(para)){
    temp.degree <- para[i,1]
    temp.span <- para[i,2]
    pre <- NULL
    for(j in 2:(length(x) - 1)) {
      y1 <- y[-j]
      x1 <- x[-j]
      model <- loess(y1 ~ x1, span = temp.span, degree = temp.degree)
      pre[j] <- predict(object = model,
                        newdata = data.frame(x1 = x[j]))
    }
    pre.all[[i]] <- pre[!is.na(pre)]
  }
  
  mse <- unlist(lapply(pre.all, function(x) {sum((y[2:(length(y)-1)] - x)^2)/length(y) }))
  result <- data.frame(para, mse, stringsAsFactors = FALSE)
  
  y.lim1 <- 0.8*min(c(y, unlist(lapply(pre.all, min))))
  y.lim2 <- 1.2*max(c(y, unlist(lapply(pre.all, max))))

  localFileName <- file.path(path, paste0("Pred_MSE_plot_", pol, ".pdf"))
  if (!file.exists(localFileName)) {
    pdf(localFileName,
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
      points(x = x[2:(length(x)-1)], y = pre.all[[k]],
             pch = 19,
             cex = 1,
             type = "o",
             col = col.list[k])
    }
    
    legend("topleft",
           legend = paste(paste("Degree", result[, 1], ", Span", result[,2]), round(result[,3], 3), sep = ": "),
           col = col.list,
           pch = 19,
           lty = 1,
           bty = "n",
           pt.cex = 1,
           lwd = 1,
           title = "Leave One Out MSE",
           cex = 1)
    dev.off()
  }
  return(result)
}