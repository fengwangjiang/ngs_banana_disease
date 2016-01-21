###############################################################################
### code chunk number : volcano plots
###############################################################################
my.volcano.plot <- function(fold.change, 
                            pvalue, 
                            lfc.threshold = 1,
                            significance = 0.001,
                            main = "Volcano plot condition day 2", 
                            is.log.fc = TRUE)
{
        if (is.log.fc)
                lfc = fold.change
        else
                lfc = log2(fold.change)
        y = -log10(pvalue)
        chosen <- (y<50)
#         lfc <- lfc[chosen]
#         y <- y[chosen]
#         pvalue <- pvalue[chosen]
        plot(lfc, y, pch = 16, col = "gray32", main = main, cex=0.45,
             xlab = "log2 fold change",
             ylab = "-log10 p-value",
             xlim = c(-max(abs(lfc)), max(abs(lfc)))
        )
        tmp = cbind(lfc, pvalue)
        high.light = abs(tmp[,1]) >= abs(lfc.threshold) &
                tmp[,2] < significance
        tmp = tmp[high.light, ]
        x = tmp[,1]
        y = -log10(tmp[,2])
        points(x = x, y = y, pch = 16, col = "red3", cex=0.45)
        abline(h = -log10(significance), col = "red", lty = 2)
        abline(v = lfc.threshold, col = "blue", lty = 2)
        abline(v = -lfc.threshold, col = "blue", lty = 2)
}
my.volcano.plot.pdf = function (fold.change, 
                                pvalue, 
                                filename,
                                lfc.threshold = 1,
                                significance = 0.001,
                                main = "Volcano plot condition day 2", 
                                is.log.fc = TRUE)
{
        # filename <- file.path(DIR.FIGURE, filename)
        pdf(file = filename)
        my.volcano.plot(fold.change, 
                        pvalue, 
                        lfc.threshold,
                        significance,
                        main, 
                        is.log.fc
        )
        tmp = dev.off()
}

###############################################################################
### code chunk number : helper functions
###############################################################################
my.volcano.plot.pdf.wrapper <- function(res, filename, main){
        lfc <- res$log2FoldChange
        pval <- res$pvalue
        my.volcano.plot.pdf(fold.change = lfc, pvalue = pval,
                            main = main, filename = filename)
}
