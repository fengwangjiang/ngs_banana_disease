KEGG.DIR <- "~/banana/musa_kegg/pathways"
files <- list.files(path = KEGG.DIR, pattern = "xml$")
files <- file.path(KEGG.DIR, files)
# pathway.id <- "01100"
# PATHWAY <- paste0("mus", pathway.id)
# PATHWAY.FILE <- paste0(PATHWAY, ".xml")
# xml.file=file.path(KEGG.DIR, PATHWAY.FILE)
for (f in files){
        node.data=node.info(f)
        if (any(is.na(unique(node.data$width))))
                message("Pathway: ", sprintf("%s", basename(f)), " has NA width.")
}
# Pathway: mus00510.xml has NA width.
# Pathway: mus00511.xml has NA width.
# Pathway: mus00514.xml has NA width.
# Pathway: mus00531.xml has NA width.
# Pathway: mus00563.xml has NA width.
# Pathway: mus00603.xml has NA width.
# Pathway: mus00604.xml has NA width.
# Pathway: mus01040.xml has NA width.
# Pathway: mus01100.xml has NA width.
# Pathway: mus01110.xml has NA width.
# Pathway: mus01200.xml has NA width.
# Pathway: mus01210.xml has NA width.
# Pathway: mus01212.xml has NA width.
# Pathway: mus01220.xml has NA width.
# Pathway: mus01230.xml has NA width.