my.go.analysis <- function(res, GO="MF", gs.pvalue=0.05, gs.size=10) {
        source("my_runGSA.R")
        # at most gs.size (10) gene sets with significant level gs.pvalue (5%)
        musa_mart_dir <- "/Users/jason/banana/banana_genome_hub/musa_mart"
        GOs <- c("MF", "BP", "CC")
        GO <- pmatch(GO, GOs)
        if(is.na(GO)) stop("invalid GO, should be MF, BP or CC")
        GO <- GOs[GO]
        filename <- paste0("gene_", GO, "_GO_term.txt")
        filename <- file.path(musa_mart_dir, filename)
        gene_MF_GO_term <- read.csv(file = filename, header = TRUE,
                                    sep = ",", quote = "\"")
        gene_MF_GO_term <- 
                dplyr::mutate(gene_MF_GO_term,
                              Gene.Uniquename=sub("P", "G", Gene.Uniquename))
#         Gene.Uniquename               Name    GO.Term             Description
#         GSMUA_Achr10G27680_001 molecular_function GO:0003677    DNA binding
        gsaRes <- my.runGSA(res = res, gene_GO = dplyr::select(gene_MF_GO_term, 
                                                    c(Gene.Uniquename, GO.Term)))
        gsaResTab <- GSAsummaryTable(gsaRes)
        # Which columns contain p-values:
#         grep("p \\(",colnames(gsaResTab),value=T)
#         grep("p \\(",colnames(gsaResTab))
#         length(gsaResTab[,10])
#         quantile(gsaResTab[,10], probs=c(seq(0.01,0.1,0.01), seq(0.2, 1, 0.1)))
        nGeneSets <- gsaRes$info$nGeneSets
        gs.size.percent <- gs.size / nGeneSets
        pvalues <- gsaResTab[, 10]
        # at most gs.size (10) gene sets with significant level gs.pvalue (5%)
        thresh <- min(gs.pvalue, quantile(pvalues, probs=gs.size.percent))
        ii <- which(pvalues < thresh)
        sig.geneset <- gsaResTab$Name[ii]
        GO_term <- unique(gene_MF_GO_term %>%
                                  dplyr::select(GO.Term, Description) %>%
                                  dplyr::filter(GO.Term %in% sig.geneset))
        GO_pvalue <- gsaResTab[ii, c(1, 10)]
        colnames(GO_pvalue) <- c("GO.Term", "pvalue")
        GO_term_pvalue <- dplyr::left_join(GO_term, GO_pvalue, by="GO.Term")
        message(sprintf("significe level: %.2f\n", thresh))
        return(GO_term_pvalue)
}
# res <- res_list_d2[[5]]
# 
# GO_term_MF <- my.go.analysis(res, GO = "MF", gs.pvalue=0.05, gs.size=10)
# GO_term_BP <- my.go.analysis(res, GO = "BP", gs.pvalue=0.05, gs.size=10)
# GO_term_CC <- my.go.analysis(res, GO = "CC", gs.pvalue=0.05, gs.size=10)
# GO_term_MF
# GO_term_BP
# GO_term_CC
# write.table(GO_term_CC, sep = "\t", row.names = FALSE)
my.go.analysis.write.table <- function(res, filename, GO="MF", ...) {
        GO_term <- my.go.analysis(res, GO = GO, ...)
        d <- dirname(filename)
        f <- basename(filename)
        f <- paste0(GO, "_", f)
        filename <- file.path(d, f)
        write.table(GO_term, file = filename, sep = "\t", row.names = FALSE)
        invisible(GO_term)
}
