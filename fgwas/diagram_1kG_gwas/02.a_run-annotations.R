
"%&%" <- function(a,b) paste0(a,b)
library("data.table")
library("dplyr")
library(GenomicRanges)

serv.dir <- "/well/got2d/jason/"
#serv.dir <- "/Users/jtorres/FUSE/"
gwas.dir <- serv.dir %&% "reference/gwas/diagram_1Kgenomes/"

mychrom <- commandArgs(trailingOnly = TRUE)[1]
if (!grepl("chr", mychrom)) {
  mychrom <- "chr" %&% mychrom
}

mypart <- commandArgs(trailingOnly = TRUE)[2]

#cred.dir <- gwas.dir %&% "DIAGRAM_T2D_Metabochip_1000G_CredibleSets_Gaulton_2015/"
#cred.files <- list.files(cred.dir)[grepl(".out",list.files(cred.dir))]
#names(cred.files) <- 1:length(cred.files)

bed.dir <- serv.dir %&% "reference/islet/"

save.dir <- serv.dir %&% "projects/t2d-integration/fgwas/diagram_1kG_gwas/intermediate_files/"
rds.dir <- serv.dir %&% "projects/t2d-integration/fgwas/diagram_1kG_gwas/rds/"

fname <- save.dir %&% mychrom %&% ".p" %&% mypart %&% ".fgwas-core.txt.gz"
core.df <- fread("cat " %&% fname %&% " | zmore")

g <- function(df){
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)>=6){
    df$strand <- gsub(pattern="[^+-]+",replacement='*', df$strand)
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}

get_overlaps <- function(a,b){
  # a and b are genomic ranges
  val <- as.integer(a %over% b)
  return(val)
}


get_annot <- function(core.df, annot, prefix, annot.input.df=TRUE){
  pb <- txtProgressBar(min = 0, max = dim(core.df)[1], initial = 0, style = 3)
  vec <- as.integer(sapply(1:dim(core.df)[1], function(i){
    setTxtProgressBar(pb, i)
    chr <- core.df$CHR[i]
    pos <- core.df$POS[i]
    snp.gr <- GRanges(chr,IRanges(pos, pos))
    if (annot.input.df==TRUE){
      annot.gr <- g(annot)
      value <- get_overlaps(snp.gr,annot.gr)
    } else{
      value <- get_overlaps(snp.gr,annot)
    }
    return(value)
  }))
  close(pb)
  saveRDS(vec,file=rds.dir%&%mychrom%&%"_"%&%mypart%&%"_"%&%prefix%&%".RDS")
  return(vec)
}

# Islet ATAC and Chromatin State Annotations

## Islet ATAC-seq

prepare_islet_atac <- function(){
  atac.df <- fread(bed.dir %&% "islet_atac_peaks.bed")
  names(atac.df) <- c("chr","start","end","id")
  atac.df$start <- atac.df$start+1
  atac.df$end <- atac.df$end+1
  return(atac.df)
}

annot_islet_atac <- function(){
  atac.df <- prepare_islet_atac()
  islet_atac <- get_annot(core.df,atac.df, "islet_atac")
  return(islet_atac)
}

islet_atac <- annot_islet_atac()


## Process Islet 15 chromatin state file


prepare_islet_chromHMM_15states <- function(){
  chmm.df <- read.table(bed.dir %&% "Pancreat_islet_15_dense.reformatted_colours.bed", sep="\t",
                        header=FALSE,stringsAsFactors = FALSE)
  chmm.df <- select(chmm.df, one_of("V1","V2","V3","V4"))
  names(chmm.df) <- c("chr","start","end","id")
  chmm.df$start <- chmm.df$start+1
  chmm.df$end <- chmm.df$end+1
  return(chmm.df)
}

annot_chrom_states <- function(){
  chmm.df <- prepare_islet_chromHMM_15states()
  for (i in 1:15){
    pre <- "islet_state" %&% i
    print(pre)
    sub.df <- filter(chmm.df,id==i)
    annot <- get_annot(core.df,sub.df, pre)
  }
}

annot_chrom_states()

# TxDb for genomic features


library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
head(seqlevels(txdb))

columns(txdb)

# Coding sequence hg 19

annot_cds <- function(){
  cds.gr <- cds(txdb)
  get_annot(core.df,cds.gr, "cds", annot.input.df = FALSE)
}

annot_transcript <- function(){
  trans.gr <- transcripts(txdb)
  get_annot(core.df,trans.gr, "transcript", annot.input.df = FALSE)
}

annot_exon <- function(){
  exon.gr <- exons(txdb)
  get_annot(core.df,exon.gr, "exon", annot.input.df = FALSE)
}

annot_intron <- function(){
  intron.gr <- unlist(intronsByTranscript(txdb))
  get_annot(core.df,intron.gr, "intron", annot.input.df = FALSE)
}

annot_5utr <- function(){
  utr5.gr <- unlist(fiveUTRsByTranscript(txdb))
  get_annot(core.df,utr5.gr, "utr_5", annot.input.df = FALSE)
}

annot_3utr <- function(){
  utr3.gr <- unlist(threeUTRsByTranscript(txdb))
  get_annot(core.df,utr3.gr, "utr_3", annot.input.df = FALSE)
}

annot_promoters <- function(){
  prom.gr <- promoters(txdb)
  get_annot(core.df,prom.gr, "promoter", annot.input.df = FALSE)
}

annot_microRNAs <- function(){
  library("mirbase.db")
  mic.gr <- microRNAs(txdb)
  get_annot(core.df,mic.gr, "micro_rna", annot.input.df = FALSE)
}

annot_dist_tss <- function(core.df,prefix){
  print("Getting Distance to Nearest TSS for each SNP")
  chr <- core.df$CHR
  pos <- core.df$POS
  snp.gr <- GRanges(chr,IRanges(pos, pos))
  trans <- transcripts(txdb)
  tss <- resize(trans, width=1, fix='start') # get TSS from transcripts
  dist <- distanceToNearest(x=snp.gr,subject=tss)
  vec <- elementMetadata(dist)$distance
  #saveRDS(vec,file=save.dir%&%prefix%&%".RDS")
  saveRDS(vec,file=rds.dir%&%mychrom%&%"_"%&%mypart%&%"_"%&%prefix%&%".RDS")
}


annot_cds()
annot_transcript()
annot_exon()
annot_intron()
annot_5utr()
annot_3utr()
annot_promoters()
##annot_microRNAs()
annot_dist_tss(core.df,"distance_tss")

## Annotation Hub

#library(AnnotationHub)
#ah = AnnotationHub()

## Epigenome RoadMap, adult pancreas and muscle DHS Narrow Peaks (mac2)

#epiFiles <- query(ah, c("GRanges","EpigenomeRoadMap","Homo sapiens"))
#df <- mcols(epiFiles)
#unique(epiFiles$species)
#unique(epiFiles$genome)
#table(epiFiles$sourcetype)
#sort(table(epiFiles$description), decreasing=TRUE)

#dhs <- query(ah , c("EpigenomeRoadMap", "GRanges", "Narrow DNasePeaks for consolidated epigenomes from EpigenomeRoadMap Project"))
#"Broad domains on enrichment for DNase-seq for consolidated epigenomes from EpigenomeRoadMap Project"
#dhs.df <- mcols(dhs)
#brain <- query(dhs,"Brain"); mcols(brain)$tags
#mus <- query(dhs,"Muscle"); mcols(mus)$tags
#panc <- query(dhs,"Pancreas"); mcols(panc)$tags

#panc_macs2.gr <- dhs[["AH30612"]]  # E098-DNase.macs2.narrowPeak.gz
#muc_macs2.gr <- dhs[["AH30627"]] #  E100-DNase.macs2.narrowPeak.gz

#annot_annhub <- function(annhub, ahid, prefix){
#  anno.gr <- annhub[[ahid]]
#  get_annot(core.df, anno.gr, prefix, annot.input.df = FALSE)
#  return("Complete")
#}

#annot_annhub(dhs,"AH30612","dhs_pancreas")
#annot_annhub(dhs,"AH30627","dhs_muscle_psoas")
