my.pdata.list <- function(dds_day){
        library(dplyr)
        library(magrittr)
        library(DESeq2)
        d <- as.data.frame(colData(dds_day))
        pdata_cell <- d %>% mutate(sampleid = rownames(d)) %>%
                dplyr::filter(condition == levels(condition)[1]) %>%
                mutate(GROUP = ifelse(cell == levels(cell)[1], 0, 1)) %>%
                dplyr::select(sampleid, GROUP)
        
        pdata_cell_in <- d %>% mutate(sampleid = rownames(d)) %>%
                dplyr::filter(condition == levels(condition)[2]) %>%
                mutate(GROUP = ifelse(cell == levels(cell)[1], 0, 1)) %>%
                dplyr::select(sampleid, GROUP)
        
        pdata_cond <- d %>% mutate(sampleid=rownames(d)) %>%
                dplyr::filter(cell==levels(cell)[1]) %>%
                mutate(GROUP=ifelse(condition==levels(condition)[1], 0, 1)) %>%
                dplyr::select(sampleid, GROUP)
        
        pdata_cond_bcl <- d %>% mutate(sampleid=rownames(d)) %>%
                dplyr::filter(cell==levels(cell)[2]) %>%
                mutate(GROUP=ifelse(condition==levels(condition)[1], 0, 1)) %>%
                dplyr::select(sampleid, GROUP)
        
        pdata_inter <- d %>% mutate(sampleid=rownames(d)) %>%
                dplyr::filter(cell==levels(cell)[1] & condition==levels(condition)[1]
                              | cell==levels(cell)[2] & condition==levels(condition)[2]) %>%
                mutate(GROUP=ifelse(cell==levels(cell)[1] & condition==levels(condition)[1], 0, 1)) %>%
                dplyr::select(sampleid, GROUP)
        # return a list of file names
        pdata_list <- list(pdata_cond, pdata_cond_bcl, pdata_cell, pdata_cell_in, pdata_inter)
#         fpdata_list <- list(fpdata_cond, fpdata_cond_bcl, fpdata_cell, fpdata_cell_in, fpdata_inter)
#         mapply(write.table, x=pdata_list, file=fpdata_list, 
#                MoreArgs = list(sep = "\t", row.names = FALSE, col.names = FALSE))
#         # return(fpdata_list)
}
# object=dds_day
# model.matrix(design(object), data = colData(object))