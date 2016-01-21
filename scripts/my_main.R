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
### code chunk number : rlog for exploratory analysis
###############################################################################
rld_file <- file.path(DIR.DATA, "rld.rds")
if (!file.exists(rld_file)){
        rld <- rlog(dds, fitType = 'local')
        saveRDS(rld, file = rld_file)
} else{
        rld <- readRDS(file = rld_file)
}
#####################################################
source("my_name_list_meta.R")

my.fname.list <- my.fname.list.meta(disease_or_drought = disease_or_drought)
my.main.list <- my.main.list.meta(disease_or_drought = disease_or_drought)
eb.dir.list <- my.eb.dir.list(disease_or_drought = disease_or_drought)
####################################################
source("my_gene_analysis_day.R")
source("my_go_analysis_day.R")
source("my_kegg_analysis_day.R")
for (day in days) {
        my.dir.create(file.path(DIR.DATA, day))
        my.gene.analysis.day(day = day)
        my.go.analysis.day(day = day)
        my.kegg.analysis.day(day = day, eb.dir.list = eb.dir.list)
}
