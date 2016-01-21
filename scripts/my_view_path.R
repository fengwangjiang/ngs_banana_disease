
my.view.path <- function (s, eset, out.prefix) 
{
        org <- substring(s, 1, 3)
        pwy.id <- sub("^[a-z]{3}", "", s)
        fc <- fData(eset)[, config.ebrowser("FC.COL")]
        gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
        gnam <- apply(fData(eset)[, gn.cols], 1, paste, collapse = ": ")
        names(fc) <- names(gnam) <- featureNames(eset)
        out.files <- paste(out.prefix, c("kview.html", "kpath.png", 
                                         "kgraph.pdf"), sep = "_")
        if (!file.exists(out.files[1])) {
                source("my_make_kpath_html.R")
                source("my_make_kgraph_html.R")
#                 debug(my.make.kpath.html)
#                 debug(my.make.kgraph.html)
                # kpath.html <- make.kpath.html(fc, pwy.id, org, out.files[2])
                kpath.html <- my.make.kpath.html(fc, pwy.id, org, out.files[2])
                # if (kpath.html==1){return(1)}
#                 kgraph.html <- make.kgraph.html(fc, gnam, pwy.id, org, 
#                                                 out.files[3])
                kgraph.html <- my.make.kgraph.html(fc, gnam, pwy.id, org, 
                                                out.files[3])
                
                cont <- EnrichmentBrowser:::make.view(kgraph.html, kpath.html, gene.html.pos = "topright")
                cat(cont, file = out.files[1])
        }
        views <- basename(out.files[1])
        return(views)
}