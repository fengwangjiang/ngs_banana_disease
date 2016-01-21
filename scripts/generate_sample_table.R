generate.sample.table <- function(disease_or_drought="disease") {
        if (disease_or_drought=="disease") {
                # generate sampleTable for banana disease experiment.
                
                SampleName <- paste0(rep(c("A", "B", "C"), each=8), rep(1:8, 3))
                condition <- rep(rep(c("Ct", "In"), each=4), 3)
                cell <- rep(rep(c("Cav", "Bcl"), each=2), 6)
                day <- rep(rep(c("D2", "D14")), 12)
                replicate <- paste0("S", rep(1:3, each=8))
                sampleTable <- data.frame(SampleName, condition, cell, day, replicate)
                rownames(sampleTable) <- SampleName
                
                write.table(sampleTable, file = EXP.DESIGN.TABLE, sep = "\t", row.names = FALSE)
        }else if (disease_or_drought=="drought"){
                # generate sampleTable for banana drought experiment.
                SampleName <- paste0(rep(c("A", "B", "C"), each=8), rep(1:8, 3))
                condition <- rep(rep(c("Wtr", "Drt"), each=4), 3)
                cell <- rep(rep(c("Cav", "Bcl"), each=2), 6)
                day <- rep(rep(c("D6", "D8")), 12)
                replicate <- paste0("S", rep(1:3, each=8))
                sampleTable <- data.frame(SampleName, condition, cell, day, replicate)
                rownames(sampleTable) <- SampleName
                
                write.table(sampleTable, file = EXP.DESIGN.TABLE, sep = "\t", row.names = FALSE)
        } else {
                stop(sprintf("disease_or_drought should be disease or drought. \n
                             Your input is %s\n", diseadisease_or_drought))
        }
}

# meta <- read.delim(file = fname)
# typeof(meta)
# class(meta)
# str(meta)