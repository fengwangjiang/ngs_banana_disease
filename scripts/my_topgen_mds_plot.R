###############################################################################
### code chunk number : top 40 rld assay for mds plot
###############################################################################
my.topgen.mds.plot <- function(res, rld, n, ...) {
        chosen <- rownames(res)[order(res$padj, -res$log2FoldChange)][1:n]
        rld <- rld[chosen, ]
        my.mds.plot(rld, ...)
}


# rld, method="euclidean", dim=2, main=NULL, filename=NULL


