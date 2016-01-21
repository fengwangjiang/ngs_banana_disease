dds_exam <- function(dds, normalized=FALSE, ...) {
        op <- par(no.readonly = TRUE)
        par(mar=c(8, 4, 2, 2))
        cnts <- counts(dds, normalized=normalized)
        cnts <- log10(cnts + 1) # avoid zeros
        boxplot(cnts, xlab="", ylab="log10 counts", las=2, ...)
        par(op)
}
# dds_exam(dds, main="original dds counts", las=2)
# dds_exam(dds_filter(dds), main="filtered dds counts", las=2)

rld_exam <- function(rld, ...) {
        op <- par(no.readonly = TRUE)
        par(mar=c(8, 4, 2, 2))
        cnts <- assay(rld)
        boxplot(cnts, xlab="", ylab="rlog of counts", las=2, ...)
        par(op)
}