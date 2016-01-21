my.gene.analysis.day <- function(day="D14"){
        ###############################################################################
        ### code chunk number : time points day2 and day14
        ###############################################################################
        my.day.subset <- function(rld_or_dds, day=levels(colData(rld_or_dds)$day)) {
                chosen <- colData(rld_or_dds)$day %in% day
                rld_or_dds_chosen <- rld_or_dds[, chosen]
                return(rld_or_dds_chosen)
        }
        
        rld_day <- my.day.subset(rld_or_dds = rld, day = day)
        dds_day <- my.day.subset(rld_or_dds = dds, day = day)
        
        # ###############################################################################
        # ### code chunk number : sample distances heatmap
        # ###############################################################################
        source("my_sampledist_pheatmap.R")
        source("mydist.R")
        for (m in METHODS) {
                my.sampledist.pheatmap(rld = rld, method = m)
                my.sampledist.pheatmap(rld = rld_day, method = m)
        }
        
        ###############################################################################
        ### code chunk number : MDS plots
        ###############################################################################
        source("my_mds_plot.R")
        for (m in METHODS) {
                for (d in c(2,3)) {
                        my.mds.plot(rld = rld, method = m, dim = d)
                        my.mds.plot(rld = rld_day, method = m, dim = d)
                }
        }
        
        ###############################################################################
        ### code chunk number : daywise data Differential Expression
        ###############################################################################
        my.day.dds.preproc <- function(dds){
                design(dds) <- ~cell + condition + cell:condition
                dds$day <- droplevels(dds$day)
                return(dds)
        }
        source("my_results.R")
        dds_day <- my.day.dds.preproc(dds_day)
        saveRDS(dds_day, file = file.path(DIR.DATA, day, "dds.rds"))
        
        res_list_day <- my.results(dds_day, alpha=0.05)
        saveRDS(res_list_day, file = file.path(DIR.DATA, day, "results_list.rds"))
        tmp_ <- sapply(res_list_day, summary)
        # list(res_cond, res_cond_bcl, res_cell, res_cell_in, res_inter)
        
        # source("my_venn.R")
        # sapply(res_list_d2, my.significant.res)
        # sapply(res_list_d14, my.significant.res)
        ###############################################################################
        ### code chunk number : day 2, 14 MA plot
        ###############################################################################
        source("my_plotMA.R")
        fname.list <- paste0("ma_", my.fname.list(day=day))
        fname.list <- file.path(DIR.FIGURE, fname.list)
        
        main_list <- paste0("MA plot ", my.main.list(day = day))
        
        tmp_ <- mapply(my.plotMA, res_list_day, fname.list, main_list)
        
        ###############################################################################
        ### code chunk number : day 2, 14 volcano plot
        ###############################################################################
        source("my_volcano.R")
        fname.list <- paste0("volcano_", my.fname.list(day=day))
        fname.list <- file.path(DIR.FIGURE, fname.list)
        
        main_list <- paste0("Volcano plot ", my.main.list(day = day))
        
        tmp_ <- mapply(my.volcano.plot.pdf.wrapper, res_list_day, 
                       fname.list, main_list)
        
        ###############################################################################
        ### code chunk number : top 40 DE genes
        ###############################################################################
        source("my_topgenes.R")
        
        fname.list <- paste0("top", NUM.TOPGENES, "genes_", my.fname.list(day=day))
        fname.list <- sub(".pdf", ".tsv", fname.list)
        fname.list <- file.path(DIR.DE, fname.list)
        
        tmp_ <- mapply(my.topgenes, res_list_day, fname.list, MoreArgs = list(n=NUM.TOPGENES))
        
        ###############################################################################
        ### code chunk number : get the top gene products
        ###############################################################################
        
        fname.list <- paste0("top", NUM.TOPGENES, "genes_function_", my.fname.list(day=day))
        fname.list <- sub(".pdf", ".tsv", fname.list)
        fname.list <- file.path(DIR.DE, fname.list)
        
        tmp_ <- mapply(my.topgeneprods, res_list_day,
                       fname.list, MoreArgs = list(n=NUM.TOPGENES))
        ###############################################################################
        ### code chunk number : top genes assay, heatmap plot
        ###############################################################################
        source("my_topgene_heatmap.R")
        method="euclidean"
        fname.list <- paste0("hm_top", NUM.TOPGENES, "genes_",  
                             method, "_", my.fname.list(day=day))
        fname.list <- file.path(DIR.FIGURE, fname.list)
        
        main_list <- paste("Top", NUM.TOPGENES, "genes heatmap\n",
                           method, my.main.list(day=day))
        
        tmp_ <- mapply(
                my.topgene.heatmap, res = res_list_day,
                filename = fname.list,
                main = main_list,
                MoreArgs = list(n = NUM.TOPGENES, method = method, rld = rld_day))
        
        ###############################################################################
        ### code chunk number : top 40 rld assay for mds plot
        ###############################################################################
        source("my_topgen_mds_plot.R")
        
        
        
        for(d in c(2, 3)){
                fname.list <- paste0("mds_top", NUM.TOPGENES, "genes_",  "dim", d, "_",
                                     method, "_", my.fname.list(day=day))
                fname.list <- file.path(DIR.FIGURE, fname.list)
                
                main_list <- paste0("Top ", NUM.TOPGENES, "genes MDS ", d, "D\n", " ", 
                                    method, " ", my.main.list(day=day))
                tmp_ <- mapply(
                        my.topgen.mds.plot, res = res_list_day,
                        filename = fname.list,
                        main = main_list,
                        MoreArgs = list(
                                dim = d, method = method,
                                rld = rld_day, n = NUM.TOPGENES))
        }
#         
#         ###############################################################################
#         ### code chunk number : Gene Set Enrichment Analysis (GO analysis)
#         ###############################################################################
#         source("my_go_analysis.R")
#         fname.list <- sub(".pdf", ".tsv", paste0("GO_Term_", my.fname.list(day=day)))
#         fname.list <- file.path(DIR.GSEA, fname.list)
#         
#         tmp_ <- mapply(my.go.analysis.write.table, res_list_day, fname.list,
#                        MoreArgs = list(GO = "MF"))
#         tmp_ <- mapply(my.go.analysis.write.table, res_list_day, fname.list,
#                        MoreArgs = list(GO = "BP"))
#         tmp_ <- mapply(my.go.analysis.write.table, res_list_day, fname.list,
#                        MoreArgs = list(GO = "CC"))
#         
}