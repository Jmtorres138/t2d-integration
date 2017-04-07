"%&%" <- function(a,b) paste0(a,b)
library("data.table")
library("dplyr")

serv.dir <- "/well/got2d/jason/"
rds.dir <- serv.dir %&% "projects/t2d-integration/fgwas/diagram_1kG_gwas/rds/"
interm.dir <- serv.dir %&% "projects/t2d-integration/fgwas/diagram_1kG_gwas/intermediate_files/"
out.dir <- serv.dir %&% "projects/t2d-integration/fgwas/diagram_1kG_gwas/fgwas_input/"



annot.vec <- c("islet_atac", "islet_state1", "islet_state2", "islet_state3",
               "islet_state4", "islet_state5", "islet_state6", "islet_state7",
               "islet_state8", "islet_state9", "islet_state10", "islet_state11",
               "islet_state12", "islet_state13", "islet_state14", "islet_state15",
               "cds", "promoter", "transcript", "exon","intron","distance_tss",
               "utr_3", "utr_5")

combine_chunk <- function(mychrom, mypart){
  fname <- interm.dir %&% "chr" %&% mychrom %&% ".p" %&% mypart %&% ".fgwas-core.txt.gz"
  core.df <- fread("cat " %&% fname %&% " | zmore")
  out.df <- core.df
  for (annot in annot.vec){
    fn <- rds.dir %&% "chr" %&% mychrom %&% "_" %&% mypart %&% "_" %&% annot %&% ".RDS"
    vec <- readRDS(fn)
    out.df <- cbind(out.df,vec)
    names(out.df)[dim(out.df)[2]] <- annot
  }
  return(out.df)
}

build_chrom <- function(mychrom){
  print("Building df for chromosome: " %&% mychrom)
  out.df <- c()
  for (p in 1:20){
    print(p)
    stack.df <- combine_chunk(mychrom,p)
    out.df <- rbind(out.df,stack.df)
  }
  return(out.df)
}

write_input_file <- function(){
  out.df <- c()
  for (c in 1:22){
    stack.df <- build_chrom(c)
    out.df <- rbind(out.df,stack.df)
  }
  out.name <- out.dir %&% "diagram_1kG_fgwas_24annot.txt.gz"
  gz <- gzfile(out.name, "w")
  write.table(out.df, gz, sep="\t", row.names=F,quote=F)
  close(gz)
}


write_input_file()
