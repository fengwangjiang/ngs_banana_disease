###########################################
disease_or_drought <- "disease"
days <- c("D2", "D14")

NUM.TOPGENES = 40
GS.SIZE <- 20
METHODS <- c("euclidean", "pearson", "spearman")
###############################################################################
### code chunk number : save figure directories
###############################################################################
RESULTS <- file.path(getwd(), "..", "Results")
DIR.FIGURE <- file.path(RESULTS, "Figures") # was savedirFig
DIR.DE <- file.path(RESULTS, "DifferentialExpression") # was savedirDE
DIR.GSEA <- file.path(RESULTS, "GSEA") # was savedirGsea

DIR.DATA <- file.path(getwd(), "..", "data")
EXP.DESIGN.TABLE <- file.path(DIR.DATA, "experiment_design_no_rownames.tsv")

my.dir.create <- function(directory){
        if (!file.exists(directory)){
                dir.create(path = directory, recursive = TRUE)
        }
}
my.dir.create(RESULTS)
my.dir.create(DIR.FIGURE)
my.dir.create(DIR.DE)
my.dir.create(DIR.GSEA)
my.dir.create(DIR.DATA)

#######################################################
MUSA_KEGG <- "/Users/jason/banana/musa_kegg"
KEGG.DIR <- "/Users/jason/banana/musa_kegg/pathways"
DIR.IDMAP <- "/Users/jason/banana/idmapping"
#######################################################
source("generate_sample_table.R")
if (!file.exists(EXP.DESIGN.TABLE)){
        generate.sample.table(disease_or_drought=disease_or_drought)
}

# ###############################################################################
# ### code chunk number : DESeq data set
# ###############################################################################
source("load_dds.R")
source("dds_filter.R")
source("dds_exam.R")

dds_file <- file.path(DIR.DATA, "dds.rds")
if (! file.exists(dds_file)) {
        source("generate_dds.R")
        dds <- generate_dds(disease_or_drought=disease_or_drought)
} else {
        dds <- load_dds()
}

# dds_exam(dds = dds)
# dds_exam(dds = dds_filter(dds))
dds <- dds_filter(dds)

###############################################################################
### code chunk number : Time course
###############################################################################
col.data <- colData(dds)
col.data$cond_cell <- factor(paste0(col.data$condition, ".", col.data$cell),
                             levels = c("Ct.Cav", "Ct.Bcl", "In.Cav", "In.Bcl"))
ddsTC <- dds
colData(ddsTC) <- col.data
ddsTC$day = relevel(ddsTC$day, "D2")
design(ddsTC) <- ~ cond_cell + day + day:cond_cell
library("BiocParallel")
register(MulticoreParam(multicoreWorkers()))
ddsTC <- DESeq(ddsTC, test="LRT", reduced=~ cond_cell + day, parallel = TRUE)
resTC <- results(ddsTC)
# results(dds, contrast=list( c("A_changed_one","A_changed_two"),
#                             c("B_changed_one","B_changed_two")), listValues=c(1/2, -1/2))
head(resTC[order(resTC$padj),],4)

resultsNames(ddsTC)
# [1] "Intercept"                  "cond_cell_Ct.Bcl_vs_Ct.Cav" "cond_cell_In.Cav_vs_Ct.Cav" "cond_cell_In.Bcl_vs_Ct.Cav"
# [5] "day_D14_vs_D2"              "cond_cellCt.Bcl.dayD14"     "cond_cellIn.Cav.dayD14"     "cond_cellIn.Bcl.dayD14"    
library(ggplot2)
fname <- file.path(DIR.FIGURE, "timecourse_onegene.pdf")
data <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup=c("day","cond_cell"), returnData=TRUE)
g <- ggplot(data, aes(x=day, y=count, color=cond_cell, group=cond_cell)) + 
        geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10() +
        ggtitle("Normalized counts for a gene with condition_cell-specific changes over time.")
ggsave(fname)


## ------------------------------------------------------------------------
betas <- coef(ddsTC)
colnames(betas)

## ----fissionheatmap------------------------------------------------------
library("pheatmap")
fname <- file.path(DIR.FIGURE, "timecourse_heatmap.pdf")
mainname <- "Heatmap of log2 fold changes for genes with smallest adjusted p value"
topGenes <- head(order(resTC$padj),30)
mat <- betas[topGenes, -c(1,2, 3, 4)]
# mat <- betas[topGenes, ]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE, filename = fname, main = mainname)
# For Ct.Cav, day14 and day2 almost no difference.
# For Ct.Bcl, day14 is a little bit upregulated.
# For In.Cav, day14 is very downregulated.
# For In.bcl, day14 is a bit downregulated.
# The same information also shows up in the time course for the most significant gene.