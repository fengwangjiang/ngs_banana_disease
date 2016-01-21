my.ea.browse <- function (res, nr.show = -1, graph.view = NULL, html.only = FALSE, day="D14") 
{
        library(EnrichmentBrowser)
        library(ReportingTools)
        library(hwriter)
        method <- ifelse(is(res$method, "character"), res$method, 
                         NA)
        eset <- res$eset
        alpha <- res$alpha
        gs <- res$gs
        out.dir <- file.path(RESULTS, "enrich_browser", day)
        # out.dir <- config.ebrowser("OUTDIR.DEFAULT")
        if (!file.exists(out.dir)) 
                dir.create(out.dir, recursive = TRUE)
        REPORTS <- "reports"
        rep.dir <- file.path(out.dir, REPORTS)
        # rep.dir <- file.path(out.dir, "reports")
        
        if (!file.exists(rep.dir)) 
                dir.create(rep.dir, recursive = TRUE)
        
        if (nr.show < 1) 
                nr.show <- res$nr.sigs
        if (nr.show > nrow(res$res.tbl)) 
                nr.show <- nrow(res$res.tbl)
        res <- res$res.tbl[seq_len(nr.show), ]
        gsc <- EnrichmentBrowser:::gs.list.2.gs.coll(gs[res[, 1]])
        res[, 1] <- sapply(res[, 1], function(s) unlist(strsplit(s, 
                                                                 "_"))[1])
        is.kegg <- is(collectionType(gsc[[1]]), "KEGGCollection")
        is.go <- is(collectionType(gsc[[1]]), "GOCollection")
        gs.title <- sapply(gsc, description)
        nr.genes <- sapply(gsc, function(g) length(geneIds(g)))
        cnames <- c(colnames(res)[1], "TITLE")
        resn <- DataFrame(res[, 1], gs.title)
        if (!("NR.GENES" %in% colnames(res))) {
                cnames <- c(cnames, "NR.GENES")
                resn <- DataFrame(resn, nr.genes)
        }
        cnames <- c(cnames, colnames(res)[2:ncol(res)])
        resn <- DataFrame(resn, res[, 2:ncol(res)])
        colnames(resn) <- cnames
        res <- resn
        im <- incidence(gsc)
        org <- organism(gsc[[1]])
        if (org == "") 
                org <- annotation(eset)
        if (!length(org)) 
                stop("Organism annotation not found!\n", "Organism under study must be annotated via annotation(eset)")
#         message("Creating gene report ...")
#         eset <- eset[colnames(im), ]
        fDat <- fData(eset)[, sapply(c("FC.COL", "ADJP.COL"), config.ebrowser)]
#         gt <- suppressMessages(gene.table(im, org, fcs = fDat))
#         gn.cols <- sapply(c("SYM.COL", "GN.COL"), config.ebrowser)
#         fData(eset)[, gn.cols] <- gt[, gn.cols]
#         gt.reps <- sapply(gsc, function(s) gene.report(s, gt, out.dir))
        link <- paste0(names(gsc), ".html")
#         res[, "NR.GENES"] <- hwrite(res[, "NR.GENES"], link = link, 
#                                     table = FALSE)
        res[, "NR.GENES"] <- hwrite(res[, "NR.GENES"],
                                    table = FALSE)
#         message("Creating set view ...")
        out.prefix <- file.path(rep.dir, names(gsc))
        names(out.prefix) <- names(gsc)
#         vcol <- sapply(gsc, function(s) view.set(eset[geneIds(s), 
#                                                       ], out.prefix[setName(s)]))
#         vcol <- hwriteImage(sub("sview.html", "volc.png", vcol), 
#                             link = vcol, table = FALSE, height = 50, width = 50, 
#                             target = "_blank")
#         res <- DataFrame(res, vcol)
#         colnames(res)[ncol(res)] <- "SET.VIEW"
        if (is.kegg) {
                message("Creating kegg view ...")
                source("my_view_path.R")
                # debug(my.view.path)
#                 vcol <- sapply(gsc, function(s) EnrichmentBrowser:::view.path(setName(s), 
#                                                           eset[geneIds(s), ], out.prefix[setName(s)]))
                vcol <- sapply(gsc, function(s) my.view.path(setName(s),
                                                             eset[geneIds(s), ],
                                                             out.prefix[setName(s)]))
                vcol <- hwriteImage(sub("kview.html", "kpath.png", vcol), 
                                    link = vcol, table = FALSE, height = 50, width = 50, 
                                    target = "_blank")
                res <- DataFrame(res, vcol)
                colnames(res)[ncol(res)] <- "PATH.VIEW"
        }
        if (!is.null(graph.view)) {
                message("Creating graph view ...")
                vcol <- sapply(gsc, function(s)
                        EnrichmentBrowser:::view.graph(eset[geneIds(s), ],
                                                       EnrichmentBrowser:::query.grn(geneIds(s), graph.view, index = FALSE), 
                                                           alpha, out.prefix[setName(s)]))
                vcol <- hwriteImage(sub("html$", "png", vcol), link = vcol, 
                                    table = FALSE, height = 50, width = 50, target = "_blank")
                res <- DataFrame(res, vcol)
                colnames(res)[ncol(res)] <- "GRAPH.VIEW"
        }
        link <- NULL
        GS.COL <- config.ebrowser("GS.COL")
        if (is.kegg) 
                link <- sapply(gsc, function(s) EnrichmentBrowser:::get.html.of.marked.pathway(setName(s), 
                                                                           geneIds(s)[fDat[geneIds(s), 2] < alpha]))
        else if (is.go) 
                link <- paste0(config.ebrowser("GO.SHOW.URL"), res[, 
                                                                   GS.COL])
        if (!is.null(link)) 
                res[, GS.COL] <- hwrite(res[, GS.COL], link = link, table = FALSE)
        htmlRep <- HTMLReport(shortName = method, title = paste(toupper(method), 
                                                                config.ebrowser("RESULT.TITLE"), sep = " - "), basePath = out.dir, 
                              reportDirectory = REPORTS)
        res <- as.data.frame(res)
        publish(res, htmlRep)
        rep <- finish(htmlRep)
        if (!html.only) 
                if (interactive()) 
                        browseURL(rep)
}