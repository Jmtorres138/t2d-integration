
# Setup

args = commandArgs(trailingOnly=TRUE)
loc.id <- args[1]

"%&%" <- function(a,b) paste0(a,b)

library("data.table")
library("dplyr")
library("ggplot2")
library("GenomicRanges")

#serv.dir <- "/Users/jtorres/FUSE/"
serv.dir <- "/well/got2d/jason/"

work.dir <- serv.dir %&% "projects/t2d-integration/fgwas/diagram_hrc/ukbb-diamante-euro-manuscript/"

#fgwas.output.dir <- work.dir %&% "fgwas_output/"
fgwas.output.dir <- work.dir %&% "conditional/fgwas_output_files/" %&% loc.id %&% "/"
pre <- fgwas.output.dir %&% "fgwas_run_loci-partition"

#region.file <- work.dir%&%"t2d-loci-regions.txt"
region.file <- work.dir %&% "region_files/" %&% "t2d-loci-regions-" %&% loc.id %&% ".txt"


get_loci_blocks <- function(save.prefix){
  loc.df <- fread(region.file)
  names(loc.df) <- c("CHR","Start","End")
  loc.gr <- GRanges(seqnames=loc.df$CHR,IRanges(loc.df$Start,loc.df$End))
  blk.df <- fread("cat " %&% pre %&% ".segbfs.gz" %&% " | zmore")
  blk.gr <- GRanges(seqnames=blk.df$chr,IRanges(blk.df$st,blk.df$sp))
  keep.gr <- blk.gr[blk.gr%over%loc.gr]
  out.df <- c()
  pb <- txtProgressBar(min=0,max=dim(loc.df)[1],style=3)
  for (i in 1:length(loc.gr)){
    setTxtProgressBar(pb,i)
    build.gr <- blk.gr[blk.gr%over%loc.gr[i]]
    chrom.vec <- as.character(seqnames(build.gr))
    start.vec <- start(build.gr)
    end.vec <- end(build.gr)
    build.df <- filter(blk.df,chr%in%chrom.vec,st%in%start.vec,sp%in%end.vec)
    SEGNUMBER <- rep(i,length(chrom.vec))
    build.df <- cbind(SEGNUMBER,build.df) %>% arrange(desc(PPA))
    out.df <- rbind(out.df,build.df[1,])
  }
  sig.blk.df <- dplyr::select(out.df,one_of("SEGNUMBER","chr","st","sp","PPA","NSNP"))

  names(sig.blk.df) <- c("SEGNUMBER","CHR","blk.start","blk.end","blk.PPA","blk.NSNP")


  write.table(x=sig.blk.df,file = save.prefix %&%".txt",
              sep="\t",quote=FALSE,row.names=FALSE)
  return(sig.blk.df)
}

save.name <- fgwas.output.dir %&% "fgwas_blk"
loc.blk.df <- get_loci_blocks(save.name)



#python subset_to_seg.py pre fgwas_blk-t2d_loci_152.txt
