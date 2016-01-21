x=res_list_d2[[5]]
x.sig=my.significant.res(x)
y=res_list_d14[[5]]
y.sig=my.significant.res(y)

x.sig.10=x.sig[1:10,]
x.sig.10
dd=x.sig.10
dd[order(dd$padj), ]

x.sig.gene <- rownames(x.sig)
y.sig.gene <- rownames(y.sig)
library(VennDiagram)
venn.diagram(x = list(day2=x.sig.gene, day14=y.sig.gene),
             filename = "/tmp/venn_inter_day2_day14.png",
             imagetype = "png",
             height = 1500, width = 2000, resolution = 500,
             main = "intersection day2 and day14",
             col="transparent", 
             cat.col=c("darkorchid1", "cornflowerblue"),
             cat.pos = 0,
             # cat.dist = -0.5,
             fill=c("darkorchid1", "cornflowerblue"),
             alpha=0.5,
             margin=0.1,
             ext.dist=0.01)
both <- intersect(x.sig.gene, y.sig.gene)


geneprod <- readRDS(file.path(DIR.DATA,"acceProtGeneProd.RDS"))[, c(3, 4)]
geneprod <- unique(geneprod)
rownames(geneprod) <- geneprod[, 1]
geneprod <- geneprod[both, ]
geneprod <- na.omit(geneprod)
# geneprod[, 2] <- paste(substr(geneprod[, 2], 1, 50), "...", sep = "")
fname <- file.path(RESULTS, "DifferentialExpression/both_Inter_D2_D14_gene_prod.tsv")
write.table(geneprod, sep = "\t", file = fname, 
            row.names = FALSE,
            col.names = c("gene", "product"))
