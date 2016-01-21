###############################################################################
### code chunk number : top genes assay, heatmap plot
###############################################################################

my.topgene.heatmap <- function(res, rld, method="euclidean", n=10, ...){
        chosen <- rownames(res)[order(res$padj, -res$log2FoldChange)][1:n]
        rld <- rld[chosen, ]
        assay <- assay(rld)
        
        meta <- as.data.frame(colData(rld))
        sample.names = with(meta, paste(condition,cell,day,replicate,sep=""))
        
        annotation_col = meta[, c("condition", "cell")]
        rownames(annotation_col) = sample.names
        
        dotargs <- list(...)
#         day <- paste(as.character(unique(meta$day)), collapse = "")
#         main <- paste("hm_top", n, "genes_\n", method, " ", day, sep = "")
#         
#         filename <- paste0("hm_top", n, "genes_", method, day, ".pdf")
#         filename <- file.path(DIR.FIGURE, filename)
#         main=NA
#         filename=NA
        mat <- assay
        arglist <- list(mat=mat,# main=main, filename=filename,
                        cluster_rows = FALSE, 
                        clustering_method = "average",
                        annotation_col = annotation_col,
                        legend=FALSE,
#                         cellheight=15,
#                         cellwidth=12,
                        height=10,
                        width=8)
        arglist <- modifyList(arglist, dotargs)
        do.call("pheatmap", arglist)
}
# sapply(res_list_d2, my.topgene.heatmap, rld = rld_day2, n = 10)