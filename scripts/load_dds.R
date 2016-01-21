load_dds <- function() {
        sampleTable <- read.delim(file = EXP.DESIGN.TABLE)
        sample.names = with(sampleTable, paste(condition,cell,day,replicate,sep=""))
        
        dds <- readRDS(file = file.path(DIR.DATA, "dds.rds"))
        colnames(dds) = sample.names
        return(dds)
}
# dds <- load_dds()
# dds
