my.make.kgraph.html <- function (fc, gname, pwy.id, org, img.file) 
{
        org <- "mus"
        width <- config.ebrowser("PLOT.WIDTH")
        height <- config.ebrowser("PLOT.HEIGHT")
        out.dir <- dirname(img.file)
        # png(img.file, width = width, height = height)
        # par(mai = rep(0, 4))
#         gr <- suppressWarnings(suppressMessages(pathview2(gene.data = fc, 
#                                                           pathway.id = pwy.id, species = org, kegg.dir = out.dir)))
        # source('my_pathview2.R')
        # debug(my.pathview2)
        gr <- suppressWarnings(suppressMessages(my.pathview2(gene.data = fc, 
                                                          pathway.id = pwy.id)))
        # dev.off()
        
        
        pv.out <- file.path(getwd(), paste0(org, pwy.id, ".kgraph.pdf"))
        file.rename(from = pv.out, to = img.file)
        kgraph.html <- sub("pdf$", "html", img.file)
        cont <- hmakeTag("html", hmakeTag("body", hwriteImage(basename(img.file))))
        cat(cont, file = kgraph.html)
        return(basename(kgraph.html))
        
        
#         kgraph.html <- sub("png$", "html", img.file)
#         if (is(gr, "graph")) {
#                 nd <- nodeRenderInfo(gr)$kegg.ids
#                 nam <- sapply(names(nd), function(n) ifelse(nd[[n]][1] %in% 
#                                                                     names(gname), gname[nd[[n]][1]], nodeRenderInfo(gr)$label[[n]]))
#                 names(nam) <- names(nd)
#                 kstr <- sapply(nd, function(n) paste(paste(org, n, sep = ":"), 
#                                                      collapse = "+"), USE.NAMES = FALSE)
#                 con <- file(kgraph.html, open = "w")
#                 refs <- paste0(config.ebrowser("KEGG.GENE.URL"), kstr)
#                 biocGraph::imageMap(gr, con = con, tags = list(HREF = refs, 
#                                                                TITLE = nam, TARGET = rep("gene", length(nd))), imgname = basename(img.file), 
#                                     width = width, height = height)
#                 close(con)
#         }
#         else cat(hmakeTag("html"), file = kgraph.html)
#         return(basename(kgraph.html))
}
my.pathview2 <- function(gene.data, pathway.id, out.suffix="kgraph"){
        # KEGG.DIR <- "~/banana/musa_kegg/pathways"
        PATHWAY <- paste0("mus", pathway.id)
        PATHWAY.FILE <- paste0(PATHWAY, ".xml")
        xml.file=file.path(KEGG.DIR, PATHWAY.FILE)
#         node.data=node.info(xml.file)
#         
#         plot.data.gene=node.map(mol.data=gene.data, node.data,
#                                 node.types="gene", node.sum = "sum")
#         cols.ts.gene=node.color(plot.data.gene, limit=1, bins=10)
#         keggview.native(plot.data.gene = plot.data.gene, node.data = node.data, same.layer = TRUE,
#                         kegg.dir = KEGG.DIR,
#                         pathway.name = PATHWAY,
#                         cols.ts.gene = cols.ts.gene, out.suffix = out.suffix )
        
        gR1 = pathview:::parseKGML2Graph2(
                xml.file, genesOnly = FALSE, expand = FALSE, split.group = FALSE
        )
        
        node.data = node.info(gR1)
        plot.data.gene=node.map(mol.data=gene.data, node.data,
                                node.types="gene", node.sum = "sum")
        plot.data.gene$width <- plot.data.gene$width*1.5 # 123456789 for gene names, previously 1.2.3.4.
        cols.ts.gene=node.color(plot.data.gene, limit=1, bins=10)
        
        keggview.graph(
                plot.data.gene = plot.data.gene, cols.ts.gene = cols.ts.gene,
                node.data = node.data, path.graph = gR1,
                pathway.name = PATHWAY, same.layer = TRUE,
                out.suffix = out.suffix
        )
        
        
}