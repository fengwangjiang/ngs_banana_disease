###############################################################################
### code chunk number : top 14 result, save
###############################################################################
my.significant.res <- function(res, lfc.thresh=1, padj.thresh=0.05){
        sig <- res$padj < padj.thresh & abs(res$log2FoldChange) > lfc.thresh
        res.sig <- res[sig, ]
        return(res.sig)
}

my.topgenes <- function(res, filename, n=40){
#         chosen <- rownames(res)[order(res$padj, -res$log2FoldChange)][1:n]
#         topres <- res[chosen, ]
        # res <- my.significant.res(res)
        topres <- res[order(res$padj, -res$log2FoldChange)[1:n], ]
        # filename <- file.path(DIR.DE, filename)
        write.table(topres[c("baseMean","log2FoldChange","padj")],
                    file = filename, sep = "\t")
        topres
}

my.topgeneprods <- function(res, filename, n=40){
        chosen <- rownames(res)[order(res$padj, -res$log2FoldChange)][1:n]
        musa_mart_dir <- "/Users/jason/banana/banana_genome_hub/musa_mart"
        gene_prod_file <- file.path(musa_mart_dir, "gene_prod.txt")
        gene_prod <- read.csv(
                file = gene_prod_file, header = TRUE,
                sep = ",", quote = "\"")
        gene_prod <- 
                dplyr::mutate(gene_prod,
                              Gene.Uniquename = sub("P", "G", Gene.Uniquename))
        rownames(gene_prod) <- gene_prod[, 1]
        gene_prod <- gene_prod[chosen, ]
        gene_prod <- na.omit(gene_prod)
        # gene_prod[, 2] <- paste(substr(gene_prod[, 2], 1, 50), "...", sep = "")
        gene_prod[, 2] <- substr(gene_prod[, 2], 1, 50)
        # gene_prod <- subset(gene_prod, gene_prod[, 1] %in% chosen)
        
        # filename <- file.path(DIR.DE, filename)
        write.table(gene_prod, file = filename, sep = "\t", row.names = FALSE)
        gene_prod
}