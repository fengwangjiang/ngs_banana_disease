my.sampledist.pheatmap <- function(rld, method="euclidean") {
        library(DESeq2)
        library(pheatmap)
        assay <- assay(rld)
        d <- mydist(assay, method)
        meta <- as.data.frame(colData(rld))
        sample.names = with(meta, paste(condition,cell,day,replicate,sep=""))
        
        # annotation_row = meta[, c("condition", "cell")]
        annotation_row = meta[, colnames(meta)[c(2,3)]]
        rownames(annotation_row) = sample.names
        
        day <- paste(as.character(unique(meta$day)), collapse = "")
        main <- paste("Sample distances\n", method, " ", day, sep = "")
        
        filename <- paste0("sample_dist_", method, day, ".pdf")
        filename <- file.path(DIR.FIGURE, filename)
        
        mat <- as.matrix(d)
        arglist <- list(mat=mat, main=main, filename=filename,
                        # cluster_rows = FALSE, 
                        clustering_method = "average",
                        annotation_row = annotation_row,
                        clustering_distance_rows = d, 
                        clustering_distance_cols = d, legend=FALSE)
        # dotnames <- modifyList(..., arglist) 
        do.call("pheatmap", arglist)
}

