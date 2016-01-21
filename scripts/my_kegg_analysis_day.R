my.kegg.analysis.day <- function(day="D14", html.only=TRUE, eb.dir.list){
        library(EnrichmentBrowser)
        source("my_musa_kegg.R")
        mus.gs <- my.load.mus.gs()
        mus.grn <- my.load.mus.grn()
        
        ###################################################
        ### mus.eset
        ###################################################
        source("load_res_list_day.R")
        res_list_day <- load.res.list.day(day=day)
        source("my_ncbigene_res_list.R")
        ncbigene_res_list <- my.ncbigene.res.list(res_list_day)
        source("my_eset.R")
        mus_eset_list <- my.make.eset.list(day=day)
        
#         eb.dir.list <- paste0(c("Cond", "Cond", "Cell", "Cell", "Inter"),
#                               c("", "_bcl", "", "_in", ""))
        source("my_ea_browse.R")
        preproc.eset <- function(mus.eset, res){
                tmp <- cbind(fData(mus.eset), 
                             res[paste0("LOC", rownames(fData(mus.eset))), c("log2FoldChange", "padj")])
                tmp <- dplyr::rename(tmp, FC=log2FoldChange, ADJ.PVAL=padj)
                tmp <- dplyr::mutate(tmp, SYMBOL=ENTREZID, GENENAME=ENTREZID)
                rownames(tmp) <- tmp$ENTREZID
                fData(mus.eset) <- tmp
                return(mus.eset)
        }
        
        for (i in seq_along(eb.dir.list)){
                mus.eset <- mus_eset_list[[i]]
                res <- ncbigene_res_list[[i]]
                mus.eset <- preproc.eset(mus.eset, res)
                sbea.res <- sbea(method="ora", eset=mus.eset, gs=mus.gs, perm=0, alpha=0.05)
                nbea.res <- nbea(method="ggea", eset=mus.eset, gs=mus.gs, grn=mus.grn, alpha = 0.2)
                res.list <- list(sbea.res, nbea.res)
                comb.res <- comb.ea.results(res.list)
                day_eb_dir <- file.path(day, eb.dir.list[[i]])
                my.ea.browse(comb.res, graph.view=mus.grn, nr.show=-1, 
                             day=day_eb_dir, html.only = html.only)
        }
}
# eb.dir.list <-
#         paste0(c("Cond", "Cond", "Cell", "Cell", "Inter"),
#                c("", "_bcl", "", "_in", ""))
# my.kegg.analysis(day="D14", html.only=TRUE, eb.dir.list = eb.dir.list)

# library(EnrichmentBrowser)
# 
# day <- "D14"
# 
# source("my_musa_kegg.R")
# mus.gs <- my.load.mus.gs()
# mus.grn <- my.load.mus.grn()
# 
# ###################################################
# ### code chunk number 29: sbea
# ###################################################
# source("load_res_list_day.R")
# res_list_day <- load.res.list.day(day=day)
# source("my_ncbigene_res_list.R")
# ncbigene_res_list <- my.ncbigene.res.list(res_list_day)
# source("my_eset.R")
# mus.eset_list <- my.make.eset.list(day=day)
# 
# 
# mus.eset <- mus.eset_list[[5]]
# 
# res <- ncbigene_res_list[[5]]
# 
# 
# preproc.eset <- function(mus.eset, res){
#         tmp <- cbind(fData(mus.eset), 
#                      res[paste0("LOC", rownames(fData(mus.eset))), c("log2FoldChange", "padj")])
#         tmp <- dplyr::rename(tmp, FC=log2FoldChange, ADJ.PVAL=padj)
#         tmp <- dplyr::mutate(tmp, SYMBOL=ENTREZID, GENENAME=ENTREZID)
#         rownames(tmp) <- tmp$ENTREZID
#         fData(mus.eset) <- tmp
#         return(mus.eset)
# }
# 
# mus.eset <- preproc.eset(mus.eset, res)
# 
# sbea.res <- sbea(method="ora", eset=mus.eset, gs=mus.gs, perm=0, alpha=0.05)
# 
# # tmp <- fData(sbea.res$eset)
# # tmp <- dplyr::mutate(tmp, SYMBOL=ENTREZID, GENENAME=ENTREZID)
# # rownames(tmp) <- tmp$ENTREZID
# # fData(sbea.res$eset) <- tmp
# 
# source("my_ea_browse.R")
# # my.ea.browse(sbea.res)
# 
# ###################################################
# ### code chunk number 36: nbea
# ###################################################
# nbea.res <- nbea(method="ggea", eset=mus.eset, gs=mus.gs, grn=mus.grn, alpha = 0.2)
# # my.ea.browse(nbea.res, graph.view = mus.grn)
# res.list <- list(sbea.res, nbea.res)
# comb.res <- comb.ea.results(res.list)
# 
# my.ea.browse(comb.res, graph.view=mus.grn, nr.show=-1, day=day)
