my.go.analysis.day <- function(day="D14") {
        
        ###############################################################################
        ### code chunk number : Gene Set Enrichment Analysis (GO analysis)
        ###############################################################################
        source("my_go_analysis.R")
        fname.list <- sub(".pdf", ".tsv", paste0("GO_Term_", my.fname.list(day=day)))
        fname.list <- file.path(DIR.GSEA, fname.list)
        
        source("load_res_list_day.R")
        res_list_day <- load.res.list.day(day=day)
        
        tmp_ <- mapply(my.go.analysis.write.table, res_list_day, fname.list,
                       MoreArgs = list(GO = "MF", gs.size=GS.SIZE))
        tmp_ <- mapply(my.go.analysis.write.table, res_list_day, fname.list,
                       MoreArgs = list(GO = "BP", gs.size=GS.SIZE))
        tmp_ <- mapply(my.go.analysis.write.table, res_list_day, fname.list,
                       MoreArgs = list(GO = "CC", gs.size=GS.SIZE))
}