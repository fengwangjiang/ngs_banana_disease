###############################################################################
### code chunk number : MA plots
###############################################################################
my.plotMA = function (res, filename, main="MA plot condition day 2")
{
        pdf(file = filename)
        DESeq2::plotMA(res, main=main)
        tmp = dev.off()
}
# 
# my.fname.list <- function(day=c("D2")){
#         fname_list <- paste0(c("Cond", "Cond", "Cell", "Cell", "Inter"), 
#                              day, c("", "_bcl", "", "_in", ""), ".pdf")
# }
# my.main.list <- function(day=c("D2")){
#         paste0(c("condition day", "condition day",
#                  "cell day", "cell day",
#                  "intersection day"), 
#                paste0(" ", day),
#                c("", " bcl", "", " inoculation", ""))
# }