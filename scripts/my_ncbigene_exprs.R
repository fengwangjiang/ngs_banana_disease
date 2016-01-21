my.ncbigene.exprs <- function(dds_day) {
        dds <- dds_day
#         source("load_dds.R")
#         source("dds_filter.R")
#         dds <- load_dds()
#         # class: DESeqDataSet
#         # dim: 37604 24
#         dds <- dds_filter(dds)
#         # class: DESeqDataSet
#         # dim: 21225 24
        # DIR.IDMAP <- "/Users/jason/banana/idmapping"
        fname <- file.path(DIR.IDMAP, "ncbiProt_ncbiGene_musaProt.rds")

        library(dplyr)
        # from GSMUA_Achr10P00010_001 to GSMUA_Achr10G00010_001
        ncbiGene.2.musaGene <- readRDS(fname)[, c(2, 3)] %>%
                mutate(musa_gene=sub("P", "G", musa_prot)) %>%
                dplyr::select(-musa_prot)
        # head(ncbiGene.2.musaGene)
        # nrow(ncbiGene.2.musaGene)
        # > nrow(ncbiGene.2.musaGene)
        # [1] 24147
        dds.df <- as.data.frame(DESeq2::counts(dds))
        dds.df <- mutate(dds.df, musa_gene=rownames(dds.df))
        # head(dds.df)
        # nrow(dds.df)
        tmp <- inner_join(ncbiGene.2.musaGene, dds.df, by="musa_gene")
        # head(tmp)
        # length(tmp)
        # nrow(tmp)
        # sum(!duplicated(tmp$NCBI_gene))
        # sum(!duplicated(tmp$musa_gene))
        tmp <- tmp[!duplicated(tmp$NCBI_gene), ]
        rownames(tmp) <- tmp$NCBI_gene
        return(tmp)
}
# source("load_dds.R")
# source("dds_filter.R")
# dds <- load_dds()
# # class: DESeqDataSet
# # dim: 37604 24
# dds <- dds_filter(dds)
# # class: DESeqDataSet
# # dim: 21225 24
#
# # dds
# # head(rownames(dds))
#
# DIR.IDMAP <- "/Users/jason/banana/idmapping/"
# fname <- file.path(DIR.IDMAP, "ncbiProt_ncbiGene_musaProt.rds")
#
# library(dplyr)
# # from GSMUA_Achr10P00010_001 to GSMUA_Achr10G00010_001
# ncbiGene.2.musaGene <- readRDS(fname)[, c(2, 3)] %>%
#         mutate(musa_gene=sub("P", "G", musa_prot)) %>%
#         dplyr::select(-musa_prot)
# # head(ncbiGene.2.musaGene)
# # nrow(ncbiGene.2.musaGene)
# # > nrow(ncbiGene.2.musaGene)
# # [1] 24147
#
#
# dds.df <- as.data.frame(DESeq2::counts(dds))
# dds.df <- mutate(dds.df, musa_gene=rownames(dds.df))
# # head(dds.df)
# # nrow(dds.df)
#
#
# tmp <- inner_join(ncbiGene.2.musaGene, dds.df, by="musa_gene")
# # head(tmp)
# # length(tmp)
# # nrow(tmp)
# # sum(!duplicated(tmp$NCBI_gene))
# # sum(!duplicated(tmp$musa_gene))
# tmp <- tmp[!duplicated(tmp$NCBI_gene), ]
