my.generate.mus.gs <- function() {
        library(EnrichmentBrowser)
#         MUSA_KEGG <- "/Users/jason/banana/musa_kegg"
        pwys <- download.kegg.pathways("mus", out.dir=MUSA_KEGG, zip=TRUE)
        pwys <- file.path(MUSA_KEGG, "mus.zip")
        mus.gs <- get.kegg.genesets(pwys = pwys)
        saveRDS(mus.gs, file.path(MUSA_KEGG,"mus_gs.rds"))
}
my.load.mus.gs <- function() {
        # MUSA_KEGG <- "/Users/jason/banana/musa_kegg"
        mus.gs <- readRDS(file.path(MUSA_KEGG, "mus_gs.rds"))
}
my.generate.mus.grn <- function(){
        # MUSA_KEGG <- "/Users/jason/banana/musa_kegg"
        pwys <- file.path(MUSA_KEGG, "mus.zip")
        mus.grn <- compile.grn.from.kegg(pwys)
        saveRDS(mus.grn, file.path(MUSA_KEGG, "mus_grn.rds"))
}
my.load.mus.grn <- function(){
        # MUSA_KEGG <- "/Users/jason/banana/musa_kegg"
        mus.grn <- readRDS(file.path(MUSA_KEGG, "mus_grn.rds"))
}
