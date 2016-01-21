###############################################################################
### code chunk number : Differential Expression, results
###############################################################################

my.results <- function(dds=dds, parallel=TRUE, ...){
        library("BiocParallel")
        register(MulticoreParam(multicoreWorkers()))
        dds <- DESeq(dds, parallel = parallel)
        resnames <- resultsNames(dds)
        # design(ddsFiltD2)
        # ~cell + condition + cell:condition
        # [1] "Intercept" "cell_Bcl_vs_Cav" "condition_In_vs_Ct"  "cellBcl.conditionIn"
        res_cell <- results(dds, name = resnames[2], parallel = parallel, ...)
        res_cond <- results(dds, name = resnames[3], parallel = parallel, ...)
        res_cond_bcl <- results(dds, contrast = list(resnames[c(3, 4)]),
                                parallel = parallel, ...)
        res_cell_in <- results(dds, contrast = list(resnames[c(2, 4)]),
                               parallel = parallel, ...)
        res_inter = results(dds, name=resnames[4], parallel = TRUE, ...)
        my.na.remove <- function(res){
                return(res[! is.na(res$padj), ])
        }
        res_cond <- my.na.remove(res_cond)
        res_cond_bcl <- my.na.remove(res_cond_bcl)
        res_cell <- my.na.remove(res_cell)
        res_cell_in <- my.na.remove(res_cell_in)
        res_inter <- my.na.remove(res_inter)
        list(res_cond, res_cond_bcl, res_cell, res_cell_in, res_inter)
}