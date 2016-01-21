my.ncbigene.res.list <- function(res_list_day) {

        # DIR.IDMAP <- "/Users/jason/banana/idmapping"
        fname <- file.path(DIR.IDMAP, "ncbiProt_ncbiGene_musaProt.rds")
        
        library(dplyr)
        # from GSMUA_Achr10P00010_001 to GSMUA_Achr10G00010_001
        ncbiGene.2.musaGene <- readRDS(fname)[, c(2, 3)] %>%
                mutate(musa_gene=sub("P", "G", musa_prot)) %>%
                dplyr::select(-musa_prot)
        ncbigene_res_list <- list()
        for (i in seq(length(res_list_day))) {
                res.df <- as.data.frame(res_list_day[[i]])
                res.df <- mutate(res.df, musa_gene=rownames(res.df))
                tmp <- inner_join(ncbiGene.2.musaGene, res.df, by="musa_gene")
                tmp <- tmp[!duplicated(tmp$NCBI_gene), ]
                rownames(tmp) <- tmp$NCBI_gene
                ncbigene_res_list[[i]] <- tmp
        }
        return(ncbigene_res_list)
}
# datadir <- "../data/D14"
# fname <- file.path(datadir, "results_list.rds") 
# res_list_day <- readRDS(fname)
# x <- my.ncbigene.res.list(res_list_day)
# length(x)
# str(x)
# head(x[[5]])
