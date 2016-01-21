library(dplyr)

# Extract protein, gene, product from NCBI banana GFF.
pgpfile <- "/Users/jason/banana/ncbi_ftp/GFF/prot_gene_prod.tsv"

# blastp NCBI protein, banana genome hub protein.
np2bpfile <- "/Users/jason/banana/banana_genome_hub/musa_blastdb/scripts/"
np2bpfile <- file.path(np2bpfile, "musa_NCBI_to_banana_blastp_filtered.tsv")

# save id mapping results, ncbi protein, ncbi gene, banana genome protein
#        NCBI_prot    NCBI_gene             musa_prot
# 1 XP_009391292.1 LOC103968336 GSMUA_Achr1P00010_001
# 2 XP_009402637.1 LOC103983387 GSMUA_Achr1P00020_001
# 3 XP_009413951.1 LOC103992244 GSMUA_Achr1P00030_001
# 4 XP_009383552.1 LOC103968337 GSMUA_Achr1P00040_001
rdsfile <- file.path(DIR.IDMAP, "ncbiProt_ncbiGene_musaProt.rds")
tsvfile <- file.path(DIR.IDMAP, "ncbiProt_ncbiGene_musaProt.tsv")
# rdsfile <- "/Users/jason/banana/idmapping/ncbiProt_ncbiGene_musaProt.rds"
# tsvfile <- "/Users/jason/banana/idmapping/ncbiProt_ncbiGene_musaProt.tsv"

prot_gene_NCBI <- read.delim(pgpfile, header = FALSE,
                             stringsAsFactors = FALSE)[, c(1,2)]
colnames(prot_gene_NCBI) <- c("NCBI_prot", "NCBI_gene")

# head(prot_gene_NCBI)
# str(prot_gene_NCBI)
# 'data.frame':	41737 obs. of  2 variables:
# $ NCBI_prot: chr  "XP_009391292.1" "XP_009402637.1" "XP_009413951.1" "XP_009383552.1" ...
# $ NCBI_gene: chr  "LOC103968336" "LOC103983387" "LOC103992244" "LOC103968337" ...

np2bp <- read.delim(np2bpfile, header = FALSE,
                    stringsAsFactors = FALSE)[, 1:3]
colnames(np2bp) <- c("NCBI_prot", "musa_prot", "identity")
np2bp <- filter(np2bp, identity > 80)

# head(np2bp)
# str(np2bp)

# 'data.frame':	36740 obs. of  3 variables:
# $ NCBI_prot: chr  "XP_009379762.1" "XP_009379764.1" "XP_009379765.1" "XP_009379766.1" ...
# $ musa_prot: chr  "GSMUA_Achr1P00360_001" "GSMUA_Achr2P00410_001" "GSMUA_Achr10P15420_001" "GSMUA_Achr10P15420_001" ...
# $ identity : num  99 100 100 100 100 ...

# summary(np2bp$identity)
# table(np2bp$identity > 80)
# sum(np2bp$identity > 80) / nrow(np2bp)
# [1] 0.8247469

# prot_gene_NCBI <- head(prot_gene_NCBI, 100)
# np2bp <- head(np2bp, 100)
np_ng_bp <- inner_join(prot_gene_NCBI, np2bp, by = "NCBI_prot") %>%
        mutate(identity=NULL)
np_ng_bp <- np_ng_bp[!duplicated(np_ng_bp$musa_prot), ]

saveRDS(np_ng_bp, file = rdsfile)
# head(np_ng_bp)
# str(np_ng_bp)
# 'data.frame':	36740 obs. of  3 variables:
# $ NCBI_prot: chr  "XP_009391292.1" "XP_009402637.1" "XP_009413951.1" "XP_009383552.1" ...
# $ NCBI_gene: chr  "LOC103968336" "LOC103983387" "LOC103992244" "LOC103968337" ...
# $ musa_prot: chr  "GSMUA_Achr1P00010_001" "GSMUA_Achr1P00020_001" "GSMUA_Achr1P00030_001" "GSMUA_Achr1P00040_001" ...

write.table(np_ng_bp, file = tsvfile, sep = "\t",
            row.names = FALSE, col.names = TRUE)
# x <- read.delim(filename, stringsAsFactors = FALSE)
# all.equal(np_ng_bp, x)
