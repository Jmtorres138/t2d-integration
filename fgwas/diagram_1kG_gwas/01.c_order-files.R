"%&%" <- function(a, b) paste0(a, b)

library("dplyr")
library("data.table")

dir <- "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_1kG_gwas/intermediate_files/"

mychrom <- commandArgs(trailingOnly = TRUE)[1]
if (!grepl("chr", mychrom)) {
    mychrom <- "chr" %&% mychrom
}

fname <- dir %&% mychrom %&% ".fgwas-core.txt.gz"

df <- fread("cat " %&% fname %&% " | zmore")
df <- arrange(df, POS)

gz <- gzfile(fname, "w")
write.table(df,file=gz,sep="\t",quote=FALSE,row.names=FALSE)
close(gz)
