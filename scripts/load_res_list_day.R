load.res.list.day <- function(day="D14") {
        res_list_day <- readRDS(file.path(DIR.DATA, day, "results_list.rds"))
        return(res_list_day)
}