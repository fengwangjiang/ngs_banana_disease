generate_dds <- function(disease_or_drought="disease") {
        # Generate summarizedExperiment se for DESeq2.
        library(Rsamtools)
        if (disease_or_drought=="disease") {
                bamdir <- "/Users/jason/banana/reads_bams/data_disease/bam"
        } else if (disease_or_drought=="drought"){
                bamdir <- "/Users/jason/banana/reads_bams/data_drought/bam"
        } else {
                stop(sprintf("disease_or_drought should be disease or drought. \n
                        Your input is %s\n", diseadisease_or_drought))
        }
        # bamdir <- "/Users/jason/banana/data_disease/bam"
        fs <- list.files(bamdir, pattern = ".sorted.bam$")
        filenames <- file.path(bamdir, fs)
        bamfiles <- BamFileList(filenames, yieldSize=2000000)
        
        ##############################################################################
        ## code chunk number : sequence features
        ##############################################################################
        library(GenomicFeatures)
        
        fname = file.path(DIR.DATA, "genes.rds")
        if (!file.exists(fname)){
                gtfdir <- "/Users/jason/banana/assembly"
                gtffile <- file.path(gtfdir, "musa.gtf")
                txdb <- makeTxDbFromGFF(gtffile, format = "gtf")
                genes <- exonsBy(txdb, by = "gene")
        } else {
                genes <- readRDS(file = fname)
        }
        ##############################################################################
        ## code chunk number : count matrix construct, extreamly time consuming
        ##############################################################################
        library(GenomicAlignments)
        
        fname = file.path(DIR.DATA, "se.rds")
        if (! file.exists(fname)) {
                se <- summarizeOverlaps(
                        features = genes, reads = bamfiles,mode = "Union",
                        singleEnd = FALSE,ignore.strand = TRUE,fragments = TRUE
                )
                saveRDS(se, file = fname)
        } else {
                se <- readRDS(file = fname)
        }
        sampleTable <- read.delim(file = EXP.DESIGN.TABLE)
        colData(se) <- DataFrame(sampleTable)
        ###############################################################################
        ### code chunk number : DESeq data set
        ###############################################################################
        library(DESeq2)
        
        fname = file.path(DIR.DATA, "dds.rds")
        if (! file.exists(fname)) {
                dds <- DESeqDataSet(se, design = ~ cell + day + condition + cell:condition)
                dds$cell = relevel(dds$cell, "Cav")
                if (disease_or_drought=="disease") {
                        dds$condition = relevel(dds$condition, "Ct")
                }
                if (disease_or_drought=="drought"){
                        dds$condition = relevel(dds$condition, "Wtr")
                }
                saveRDS(dds, file = fname)
        } else {
                dds <- readRDS(file = fname)
        }
        sample.names = with(sampleTable, paste(condition,cell,day,replicate,sep=""))
        colnames(dds) = sample.names
        return(dds)
}
