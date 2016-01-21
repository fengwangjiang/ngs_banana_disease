my.fname.list.meta <- function(disease_or_drought) {
        function(day) {
                if (disease_or_drought == "disease") {
                        fname_list <- paste0(
                                c(
                                        "Cond", "Cond", "Cell", "Cell", "Inter"
                                ),
                                day, c("", "_bcl", "", "_in", ""), ".pdf"
                        )
                } else if (disease_or_drought == "drought") {
                        fname_list <- paste0(
                                c(
                                        "Cond", "Cond", "Cell", "Cell", "Inter"
                                ),
                                day, c("", "_bcl", "", "_drt", ""), ".pdf"
                        )
                } else {
                        stop(
                                sprintf(
                                        "disease_or_drought should be disease or drought. \n
                                        Your input is %s\n", diseadisease_or_drought
                                )
                                )
                }
                return(fname_list)
                }
}

my.main.list.meta <- function(disease_or_drought) {
        function(day) {
                if (disease_or_drought == "disease") {
                        main_list <- paste0(
                                c(
                                        "condition day", "condition day",
                                        "cell day", "cell day",
                                        "intersection day"
                                ),
                                paste0(" ", day),
                                c("", " bcl", "", " inoculation", "")
                        )
                } else if (disease_or_drought == "drought") {
                        main_list <- paste0(
                                c(
                                        "condition day", "condition day",
                                        "cell day", "cell day",
                                        "intersection day"
                                ),
                                paste0(" ", day),
                                c("", " bcl", "", " drought", "")
                        )
                } else {
                        stop(
                                sprintf(
                                        "disease_or_drought should be disease or drought. \n
                                        Your input is %s\n", diseadisease_or_drought
                                )
                                )
                }
                return(main_list)
                }
}

my.eb.dir.list <- function(disease_or_drought){
        if (disease_or_drought == "disease") {
                eb.dir.list <-
                        paste0(c("Cond", "Cond", "Cell", "Cell", "Inter"),
                               c("", "_bcl", "", "_in", ""))

        } else if (disease_or_drought == "drought") {
                eb.dir.list <-
                        paste0(c("Cond", "Cond", "Cell", "Cell", "Inter"),
                               c("", "_bcl", "", "_drt", ""))
        } else {
                stop(
                        sprintf(
                                "disease_or_drought should be disease or drought. \n
                                        Your input is %s\n", diseadisease_or_drought
                        )
                )
        }
        return(eb.dir.list)
}
