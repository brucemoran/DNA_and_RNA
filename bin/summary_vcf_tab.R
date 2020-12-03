#! /usr/bin/env R

library(tidyverse)
options(scipen=999)
##make summary of VCF table from parse_merged_vcf.pl
argsIn <- commandArgs(trailingOnly = TRUE)
outpath <- gsub(".tab", ".summary.tab", argsIn[1])
outRData <- gsub(".tab", ".summary.RData", argsIn[1])

if(!exists(outRData)){
  t1 <- read_tsv(argsIn[1])
  if(length(argsIn)>1){
    tot_bp <- as.numeric(argsIn[2])
  }

  ##parse per sample for variants (all)
  samp1 <- c(grep("HGVSp", colnames(t1)))+1
  sampn <- c(length(colnames(t1)))-1
  col_samp <- samp1:sampn
  per_samp <- lapply(col_samp, function(f){
    print(colnames(t1)[f])
    t2 <- t1 %>% dplyr::select(all_of(f)) %>% unlist()
    t3 <- table(t2==".")["FALSE"]
    names(t3) <- NULL
    return(t3)
  })
  names(per_samp) <- colnames(t1)[col_samp]
  t1_cons <- t1 %>% dplyr::select(Consequence)  %>% unlist() %>% table() %>% sort()
  t1_popf <- t1 %>% dplyr::select(pop_freq)  %>% unlist() %>% table() %>% sort()

  ##parse per sample for variants (MODERATE, HIGH)
  t1_mh <- t1 %>% dplyr::filter(IMPACT %in% c("MODERATE","HIGH"))
  per_samp_mh <- lapply(col_samp, function(f){
    print(colnames(t1_mh)[f])
    t2 <- t1_mh %>% dplyr::select(all_of(f)) %>% unlist()
    t3 <- table(t2==".")["FALSE"]
    names(t3) <- NULL
    return(t3)
  })
  names(per_samp_mh) <- colnames(t1_mh)[col_samp]
  t1_mh_cons <- t1_mh %>% dplyr::select(Consequence)  %>% unlist() %>% table() %>% sort()
  t1_mh_popf <- t1_mh %>% dplyr::select(pop_freq)  %>% unlist() %>% table() %>% sort()

  ##if size of region variant-called supplied
  if(exists("tot_bp")){
    tot_mb <- tot_bp/1000000
    per_samp_tmb <- lapply(per_samp,function(f){
      f/tot_mb
    })
    per_samp_mh_tmb <- lapply(per_samp_mh,function(f){
      f/tot_mb
    })
  }
  save.image(file=outRData)
}

if(exists(outRData)){
  load(outRData)
}

##report
ato <- as.data.frame(cbind(unlist(per_samp), unlist(per_samp_mh), unlist(per_samp_tmb), unlist(per_samp_mh_tmb)))
colnames(ato) <- c("total_var", "HM_var", "total_TMB", "HM_TMB")
atot <- as_tibble(ato, rownames="sample")
write_tsv(atot, path=outpath)
