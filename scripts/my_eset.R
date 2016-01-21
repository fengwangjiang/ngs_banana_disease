# construct musa expression set. 
# we need expression data, phenotype data, and feature data.

my.make.eset.list <- function(day="D14"){
        datadir <- file.path(DIR.DATA, day)
        fname <- file.path(datadir, "dds.rds")
        dds_day <- readRDS(file = fname)
        
        fname <- file.path(datadir, "results_list.rds") 
        res_list_day <- readRDS(fname)
        
        file_list <- list(
                cond = "cond.tsv", cond_bcl = "cond_bcl.tsv",
                cell = "cell.tsv", cell_in = "cell_in.tsv",
                inter = "inter.tsv")
        
        fdata_file_list <- lapply(file_list, function(x) {
                paste0("fdata_", x)
        })
        fdata_file_list <- lapply(fdata_file_list, function(file){file.path(datadir, file)})
        
        pdata_file_list <- lapply(file_list, function(x) {
                paste0("pdata_", x)
        })
        pdata_file_list <- lapply(pdata_file_list, function(file){file.path(datadir, file)})
        
        exprs_file_list <- lapply(file_list, function(x) {
                paste0("expr_", x)
        })
        exprs_file_list <- lapply(exprs_file_list, function(file){file.path(datadir, file)})
        library(dplyr)
        source("my_ncbigene_exprs.R")
        ncbigene_exprs <- my.ncbigene.exprs(dds_day = dds_day)
        # rownames(ncbigene_exprs) <- ncbigene_exprs$NCBI_gene
        
        ncbigene_musagene <- ncbigene_exprs[, c(1, 2)]
        
#         fdata_list <- sapply(res_list_day, function(res){
#                 res.df <- as.data.frame(res)
#                 res.df <- mutate(res.df, musa_gene=rownames(res.df))
#                 tmp <- inner_join(ncbigene_musagene, res.df, by="musa_gene")
#                 tmp <- tmp[!duplicated(tmp$NCBI_gene), ]
#                 fdata <- tmp %>% dplyr::select(NCBI_gene)
#                 return(fdata)
#         })
        source("my_ncbigene_res_list.R")
        ncbigene_res_list <- my.ncbigene.res.list(res_list_day)
        
        fdata_list <- sapply(ncbigene_res_list, function(x){dplyr::select(x, NCBI_gene)})
        
        source("my_pdata_list.R")
        pdata_list <- my.pdata.list(dds_day = dds_day)
        
        plen <- length(pdata_list)
        flen <- length(fdata_list)
        stopifnot(plen == flen)
        
        exprs_list <- list()
        for (i in seq(1:plen)){
                fdata <- fdata_list[[i]]
                pdata <- pdata_list[[i]]
                exprs_list[[i]] <- ncbigene_exprs[fdata, pdata[, 1]]
        }
        # remove "LOC" from LOC103983387 to get the ENTREZ ID, which are numbers
        fdata_list <- lapply(fdata_list, function(x){gsub("LOC", "", x)})
    
        mapply(write.table, x=fdata_list, file=fdata_file_list,
               MoreArgs = list(sep = "\t", row.names = FALSE, col.names = FALSE))
        
        mapply(write.table, x=pdata_list, file=pdata_file_list,
               MoreArgs = list(sep = "\t", row.names = FALSE, col.names = FALSE))
        
        mapply(write.table, x=exprs_list, file=exprs_file_list,
               MoreArgs = list(sep = "\t", row.names = FALSE, col.names = FALSE))
        
        eset_list <- mapply(EnrichmentBrowser::read.eset, 
                            exprs.file=exprs_file_list,
                            pdat.file=pdata_file_list,
                            fdat.file=fdata_file_list)
}