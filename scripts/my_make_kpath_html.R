my.make.kpath.html <- function (fc, pwy.id, org, img.file) 
{
        org <- "mus"
        out.dir <- dirname(img.file)
#         pathview(gene.data = fc, pathway.id = pwy.id, species = org, 
#                  kegg.dir = out.dir, out.suffix = "kpath")
        my.pathview(gene.data = fc, pathway.id = pwy.id, out.suffix = "kpath")
        pv.out <- file.path(getwd(), paste0(org, pwy.id, ".kpath.png"))
        file.rename(from = pv.out, to = img.file)
        kpath.html <- sub("png$", "html", img.file)
        cont <- hmakeTag("html", hmakeTag("body", hwriteImage(basename(img.file))))
        cat(cont, file = kpath.html)
        return(basename(kpath.html))
}

my.pathview <- function(gene.data, pathway.id, out.suffix="kpath"){
        # KEGG.DIR <- "~/banana/musa_kegg/pathways"
        PATHWAY <- paste0("mus", pathway.id)
        PATHWAY.FILE <- paste0(PATHWAY, ".xml")
        xml.file=file.path(KEGG.DIR, PATHWAY.FILE)
        node.data=node.info(xml.file)
        
        if (any(is.na(node.data$width)))
        {
                png.file <- sub("xml$", "png", xml.file)
                pv.out <- file.path(getwd(), paste0("mus", pathway.id, ".kpath.png"))
                file.copy(from = png.file, to = pv.out)
                return(0)
                
        }
                
        plot.data.gene=node.map(mol.data=gene.data, node.data,
                                node.types="gene", node.sum = "sum")
        cols.ts.gene=node.color(plot.data.gene, limit=1, bins=10)
        keggview.native(plot.data.gene = plot.data.gene, node.data = node.data, same.layer = TRUE,
                        kegg.dir = KEGG.DIR,
                        pathway.name = PATHWAY,
                        cols.ts.gene = cols.ts.gene, out.suffix = out.suffix )
        # gR1=pathview:::parseKGML2Graph2(xml.file, genesOnly=FALSE, expand=FALSE, split.group=FALSE)
#         gR1=pathview:::parseKGML2Graph2(xml.file, genesOnly=TRUE, expand=FALSE, split.group=FALSE)
#         node.data2=node.info(gR1)
#         
#         keggview.graph(plot.data.gene = plot.data.gene, cols.ts.gene = cols.ts.gene, 
#                        node.data = node.data2, path.graph = gR1, 
#                        pathway.name = PATHWAY, same.layer = TRUE,
#                        out.suffix = out.suffix)
        
        
}