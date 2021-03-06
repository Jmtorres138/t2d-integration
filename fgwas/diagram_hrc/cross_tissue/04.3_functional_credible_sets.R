
#Constructing functional (fgwas) credible sets from the fgwas model with the best cross validated likelihood

# Setup

"%&%" <- function(a,b) paste0(a,b)
library("data.table")
library("dplyr")
library("ggplot2")
library(GenomicRanges)

serv.dir <- "/well/got2d/jason/"

work.dir <- serv.dir %&% "projects/t2d-integration/fgwas/diagram_hrc/cross_tissue/"

cred.set.dir <- serv.dir %&% "projects/t2d-integration/fgwas/diagram_hrc/cross_tissue/credible_sets/"


get_cred <- function(dframe,cname,prob=0.99){ # 99% functional credible sets
  index <- match(cname,names(dframe))
  vec <- sort(dframe[,index],decreasing=TRUE)
  count=0
  sum=0
  for (v in vec){
    count <- count + 1
    sum <- sum + v
    if (sum >= prob){
      break
    }
  }
  return(count)
}


# Determine functional credible sets (95% and 99% Credible Sets)

get_credsets <- function(probthresh,segsnps.df){
  seg.vec <- sort(unique(segsnps.df$SEGNUMBER))
  out.df <- c()
  pb <- txtProgressBar(min=0,max=length(seg.vec),style = 3)
  for (i  in 1:length(seg.vec)){
    setTxtProgressBar(pb,i)
    seg <- seg.vec[i]
    temp.df <- filter(segsnps.df,SEGNUMBER==seg) %>% arrange(desc(PPA))
    cumppa <- sum(temp.df$PPA)
    temp.df <- temp.df[1:get_cred(temp.df,"PPA",prob=probthresh*cumppa),]
    temp.df <- dplyr::select(temp.df,-pi,-chunk,-pseudologPO,-pseudoPPA,-V)
    names(temp.df)[c(2,3,4)] <- c("SNPID","CHR","POS")
    temp.df$PPA <- temp.df$PPA/cumppa # rescale ppa to reflect proportion of cummulative sum in block
    out.df <- rbind(out.df,temp.df)
  }
  return(out.df)
}



# Get nearest genes


library(GenomicFeatures)
library(org.Hs.eg.db)
library(annotate)


annot_refGene <- function(cred.df,segsnps.df){
  seg.vec <- sort(unique(segsnps.df$SEGNUMBER))
  snp.gr <- GRanges(cred.df$CHR,IRanges(cred.df$POS, cred.df$POS))

  hg19.refseq.db <- makeTxDbFromUCSC(genome="hg19", table="refGene")
  refseq.genes<- genes(hg19.refseq.db)
  all.geneids <- elementMetadata(refseq.genes)$gene_id
  all.genesymbols <- getSYMBOL(all.geneids, data='org.Hs.eg')
  ref.df <- data.frame(gene.id=all.geneids,symbol=all.genesymbols,
                       stringsAsFactors = FALSE)
  elementMetadata(refseq.genes)$gene_id <- all.genesymbols
  df <- as.data.frame(refseq.genes)
  sub.gr <- GRanges(seqnames=df$seqnames,
                    IRanges(start=df$start,end=df$end),
                    strand=rep("*",dim(df)[1]),
                    name=df$gene_id)
  nearestGenes <- nearest(snp.gr,sub.gr)
  res <- df$gene_id[nearestGenes]
  dist <- distance(snp.gr, refseq.genes[nearestGenes])
  symbol <- res
  cred.df <- cbind(symbol,cred.df)


  # Sync symbol names for each SEGNUMBER by MOST COMMON gene
  sync.df <- c()
  pb <- txtProgressBar(min=0,max=length(seg.vec),style=3)
  for (i in 1:length(seg.vec)){
    setTxtProgressBar(pb,i)
    seg <- seg.vec[i]
    temp.df <- filter(cred.df,SEGNUMBER==seg)
    top <- names(sort(table(temp.df$symbol),decreasing=TRUE))[1]
    temp.df$symbol <- rep(top,dim(temp.df)[1])
    sync.df <- rbind(sync.df,temp.df)
  }
  cred.df <- sync.df
  return(cred.df)
}

# Run

run_process <- function(round, tissue, ppthr="0.9"){
  segsnps.df <- fread("cat " %&% cred.set.dir %&%
                        "Round"%&%round%&%"."%&%tissue%&%".loci_block_snps.bfs.txt.gz"
                      %&% " | zmore")

  cred95.df <- get_credsets(0.95,segsnps.df)
  cred99.df <- get_credsets(0.99,segsnps.df)
  cred95.df <- annot_refGene(cred95.df,segsnps.df)
  cred99.df <- annot_refGene(cred99.df,segsnps.df)
  write.table(x=cred95.df,file=cred.set.dir%&%"Round"%&%round%&%"."%&%
                tissue%&%".fgwas_credsets_95.txt",sep="\t",
              quote=FALSE,row.names=F)
  write.table(x=cred99.df,file=cred.set.dir%&%"Round"%&%round%&%"."%&%
                tissue%&%".fgwas_credsets_99.txt",sep="\t",
              quote=FALSE,row.names=F)
}


run_jobs <- function(round,tiss.vec,ppthr="0.9"){
  print("Round" %&% round)
  for (tiss in tiss.vec){
    run_process(round,tiss,ppthr)
  }
}


run_jobs("1",tiss.vec=c("adipose","muscle","liver"))#"islet",
run_jobs("2",tiss.vec=c("adipose","muscle","liver"))
