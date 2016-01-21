###############################################################################
### code chunk number : gene sets enrichment analysis
###############################################################################
my.runGSA <- function(res, gene_GO, ...) {
        library(dplyr)
        library(snowfall)
        library(piano)
        gsc <- loadGSC(gene_GO)
        res <- as.data.frame(res)
        argslist <- list(geneLevelStats=dplyr::select(res, pvalue),
                         directions=dplyr::select(res, log2FoldChange),
                         geneSetStat="mean",
                         gsc=gsc,
                         gsSizeLim=c(5, 300),
                         signifMethod="geneSampling",
                         nPerm=5000,
                         ncpus=5,
                         verbose=TRUE)
        dotslist <- list(...)
        argslist <- modifyList(argslist, dotslist)
        gsaRes <- do.call(runGSA, argslist)
        return(gsaRes)
}
# quantile(sapply(gsc$gsc, length), probs=c(seq(0.1,0.8,0.1), seq(0.9, 1, 0.01)))