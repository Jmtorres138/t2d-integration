# Usage: 
# Rscript --vanilla 06.0_conditional-create-bed.R number_of_conditional_list_file 
# Example:
# Rscript --vanilla 06.0_conditional-create-bed.R 1

"%&%" <- function(a,b) paste0(a,b)
library("data.table")
library("dplyr")
library("ggplot2")
library("GenomicRanges")
#serv.dir <- "/Users/jtorres/FUSE/"
serv.dir <- "/well/got2d/jason/"
work.dir <- serv.dir %&% "projects/t2d-integration/fgwas/diagram_hrc/ukbb-diamante-euro-manuscript/"
fgwas.output.dir <- work.dir %&% "fgwas_output/"
cred.set.dir <- serv.dir %&% "projects/t2d-integration/fgwas/diagram_hrc/ukbb-diamante-euro-manuscript/credible_sets/"
region.dir <- work.dir %&% "region_files/"
part.dir <- work.dir %&% "partition_bedfiles/"
if (!(dir.exists(region.dir))){
  dir.create(region.dir)
}
if (!(dir.exists(part.dir))){
  dir.create(part.dir)
}
size.df <- fread(work.dir%&%"hg19.chrom.sizes")
lim.df <- fread(work.dir %&% "chromosome_limits.txt")
gen.gwas.dir <- serv.dir %&% "reference/gwas/diamante-ukbb_hrc/conditioned/"
args = commandArgs(trailingOnly=TRUE)
number <- args[1] 
print(number)
#number <- 1
loci.list.file <- gen.gwas.dir %&% "list.for.credible.sets.cond" %&% number %&% ".txt"


# Create region file 

make_loc_df <- function(){
  list.df <- fread(loci.list.file)
  names(list.df) <- c("CHR","POS","LOCUS")
  list.df$CHR <- "chr"%&%list.df$CHR
  list.df$win.start <- as.integer(list.df$POS - 500000)
  list.df$win.end <- as.integer(list.df$POS + 500000)
  for (i in 1:dim(list.df)[1]){
    c <- list.df$CHR[i]
    s <- list.df$win.start[i]
    e <- list.df$win.end[i]
    if (s < 1){
      list.df$win.start[i] <- 1 
    }
    cap <- filter(lim.df,V1==c)$V3
    if (e > cap){
      list.df$win.end[i] <- cap
    }
  }
  list.df$win.start <- as.integer(list.df$win.start)
  return(list.df)
}

process_loci_gr <- function(){
  loc.df <- make_loc_df()
  loci.gr <- GRanges(seqnames=loc.df$CHR,IRanges(loc.df$win.start,loc.df$win.end))
  out.gr <- GRanges()
  pb <- txtProgressBar(min=0,max=length(loci.gr),style=3)
  for (i in 1:length(loci.gr)){
    setTxtProgressBar(pb,i)
    gr <- loci.gr[i] 
    check <- loci.gr[loci.gr %over% gr]
    chromo <- as.character(seqnames(check))[1]
    if (length(check)==2){ # None are greater than 2 
      chunk1.gr <- GRanges(seqnames=chromo,IRanges(start(check)[1],start(check)[2]))
      chunk2.gr <- GRanges(seqnames=chromo,IRanges(end(check)[1],end(check)[2]))
      out.gr <- append(out.gr,chunk1.gr)
      out.gr <- append(out.gr,chunk2.gr)
    } else{
      out.gr <- append(out.gr, gr)
    }
  }
  return(reduce(out.gr))
}

loc.df <- make_loc_df()

loc.gr <- process_loci_gr()
loc.df <- data.frame("CHR"=seqnames(loc.gr),
                     "Start"=start(loc.gr),
                     "End"=end(loc.gr))
loc.df$CHR <- as.character(loc.df$CHR)

write.table(loc.df,file=region.dir%&%"t2d-loci-regions-cond" %&% number %&%".txt",sep="\t",quote=FALSE,
            row.names=F,col.names=F)


# Create Bed file 

process_loci_gr <- function(){
  loci.gr <- GRanges(seqnames=loc.df$CHR,IRanges(loc.df$Start,loc.df$End))
  out.gr <- GRanges()
  pb <- txtProgressBar(min=0,max=length(loci.gr),style=3)
  for (i in 1:length(loci.gr)){
    setTxtProgressBar(pb,i)
    gr <- loci.gr[i] 
    check <- loci.gr[loci.gr %over% gr]
    chromo <- as.character(seqnames(check))[1]
    if (length(check)==2){ # None are greater than 2 
      chunk1.gr <- GRanges(seqnames=chromo,IRanges(start(check)[1],start(check)[2]))
      chunk2.gr <- GRanges(seqnames=chromo,IRanges(end(check)[1],end(check)[2]))
      out.gr <- append(out.gr,chunk1.gr)
      out.gr <- append(out.gr,chunk2.gr)
    } else{
      out.gr <- append(out.gr, gr)
    }
  }
  return(reduce(out.gr))
}


divide_seq <- function(chromo,start,end,width=1e6){
  vec <- seq(start,end,width)
  vec <- unique(append(vec,end))
  out.df <- c()
  for (i in 1:length(vec)){
    if (i == length(vec)){
      pass <- TRUE 
    } else if (i == (length(vec)-1)){
      build.df <- data.frame("CHR"=chromo,"Start"=vec[i],"End"=(vec[i+1]))
      out.df <- rbind(out.df,build.df)
    } else{
      #build.df <- data.frame("CHR"=chromo,"Start"=vec[i],"End"=(vec[i+1]-1))
      build.df <- data.frame("CHR"=chromo,"Start"=vec[i],"End"=(vec[i+1]))
      
      out.df <- rbind(out.df,build.df)
    }
  }
  return(out.df)
}

get_bed_df <- function(){
  loci.gr <- process_loci_gr()
  #dis.df <- data.frame("CHR"=seqnames(disjoin.gr),
  #                     "Start"=start(disjoin.gr),
  #                     "End"=end(disjoin.gr))
  lim.gr <- GRanges(seqnames=lim.df$V1,IRanges(lim.df$V2,lim.df$V3))
  
  tot.gr <- sort(sortSeqlevels(append(lim.gr,loci.gr)))
  disjoin.gr <- disjoin(tot.gr)
  dis.df <- data.frame("CHR"=seqnames(disjoin.gr),
                       "Start"=start(disjoin.gr),
                       "End"=end(disjoin.gr))
  dis.df$CHR <- as.character(dis.df$CHR)
  dis.df$chrom <- as.integer(gsub("chr","",dis.df$CHR))
  dis.df <- arrange(dis.df,chrom)
  dis.df$range <- dis.df$End - dis.df$Start
  out.df <- c()
  for (i in 1:dim(dis.df)[1]){
    if (dis.df$range<=1e6){
      build.df <- data.frame("CHR"=dis.df$CHR[i],
                             "Start"=dis.df$Start[i],
                             "End"=dis.df$End[i])
    } else{
      build.df <- divide_seq(dis.df$CHR[i],dis.df$Start[i],dis.df$End[i])
    }
    out.df <- rbind(out.df,build.df)
  }
  out.df$Start <- out.df$Start - 1
  out.df$End <- out.df$End - 1
  
  return(out.df)
}


bed.df <- get_bed_df()
bed.df$CHR <- as.character(bed.df$CHR)
write.table(bed.df,file=part.dir%&%"loci-partition-cond" %&% number %&% ".bed",sep="\t",quote=FALSE,
            row.names=F,col.names=F)
