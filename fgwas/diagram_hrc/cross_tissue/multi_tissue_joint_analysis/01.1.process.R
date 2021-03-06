
library("data.table")
library("dplyr")
library(GenomicRanges)

save.dir <- "/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/cross_tissue/multi_tissue_joint_analysis/"

command <- paste0("cat ","/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/cross_tissue/multi_tissue_joint_analysis/fgwas_input/intermediate.fgwas.gz"," | zmore")
df <- fread(command)
print("Sorting file...")
chrom <- gsub("chr","",df$CHR)
df <- cbind(chrom,df)
df <- arrange(df,chrom,POS)
df <- dplyr::select(df,-chrom)
print(str(df))

print("Removing duplicates...")
df <- df[!duplicated(df$SNPID),]
print(str(df))

# TxDb for genomic features
library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

print("Getting Distance to Nearest TSS for each SNP")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
trans <- transcripts(txdb)
snp.gr <- GRanges(df$CHR,IRanges(df$POS0, df$POS))
tss <- resize(trans, width=1, fix='start') # get TSS from transcripts
dist <- distanceToNearest(x=snp.gr,subject=tss)
distance_tss <- elementMetadata(dist)$distance
df$distance_tss <- distance_tss

write.table(df,"/well/got2d/jason/projects/t2d-integration/fgwas/diagram_hrc/cross_tissue/multi_tissue_joint_analysis/fgwas_input/ukbb_diamante-euro.fgwas",sep="	",row.names=F,quote=F)

    