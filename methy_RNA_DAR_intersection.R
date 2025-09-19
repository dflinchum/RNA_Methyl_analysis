library(openxlsx)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(VennDiagram)
library(RColorBrewer)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicRanges)
library(ChIPpeakAnno)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene











#methylation data

setwd("C:/Users/Flinchum/Documents/Tulane Docs/Lee Lab/Results/Methyl Analysis/Dean")
## Differentially methylated REGIONS
BER_minusDox_hyper <- read.csv("BER_DMR_peakAnno_Dox_hyper.csv")
BER_plusDox_hyper <- read.csv("BER_DMR_peakAnno_pDox_hyper.csv")

JN_minusDox_hyper <- read.csv("JN_DMR_peakAnno_Dox_hyper.csv")
JN_plusDox_hyper <- read.csv("JN_DMR_peakAnno_pDox_hyper.csv")

SK2_minusDox_hyper <- read.csv("SK2_DMR_peakAnno_Dox_hyper.csv")
SK2_plusDox_hyper <- read.csv("SK2_DMR_peakAnno_pDox_hyper.csv")

Tumor2911_hyper <- read.csv("DMR_peakAnno_primary_hyper.csv")
Tumor2912_hyper <- read.csv("DMR_peakAnno_recurrent_hyper.csv")



##combining into one df with all methyl values; worked well. combined all the rows, not joined
BER_hyper_merged <- bind_rows(BER_minusDox_hyper %>% dplyr::select(seqnames, start, end, diff.Methy, strand),
                              BER_plusDox_hyper %>% dplyr::select(seqnames, start, end, diff.Methy, strand))

JN_hyper_merged <- bind_rows(JN_minusDox_hyper %>% dplyr::select(seqnames, start, end, diff.Methy, strand),
                             JN_plusDox_hyper %>% dplyr::select(seqnames, start, end, diff.Methy, strand))

SK2_hyper_merged <- bind_rows(SK2_minusDox_hyper %>% dplyr::select(seqnames, start, end, diff.Methy, strand),
                              SK2_plusDox_hyper %>% dplyr::select(seqnames, start, end, diff.Methy, strand))

Tumor291_hyper_merged <- bind_rows(Tumor2911_hyper %>% dplyr::select(seqnames, start, end, diff.Methy, strand),
                                   Tumor2912_hyper %>% dplyr::select(seqnames, start, end, diff.Methy, strand))


BER_hyper_merged <- makeGRangesFromDataFrame(BER_hyper_merged, keep.extra.columns = TRUE)
JN_hyper_merged <- makeGRangesFromDataFrame(JN_hyper_merged, keep.extra.columns = TRUE)
SK2_hyper_merged <- makeGRangesFromDataFrame(SK2_hyper_merged, keep.extra.columns = TRUE)
Tumor291_hyper_merged <- makeGRangesFromDataFrame(Tumor291_hyper_merged, keep.extra.columns = TRUE)




#WT1 binding sites
setwd("D:/6_Data Files")
WT1_enriched <- read.table("WT1_enriched_regions.bed", header = FALSE, sep = "\t")

WT1_enriched <- WT1_enriched %>% rename("chr"="V1", "start"="V2", "end"="V3", "length"="V4", "abs_summit"="V5", "pileup"="V6", "-log10pvalue"="V7", "fold_enrichment"="V8", "-log10qvalue"="V9", "name"="V10")
WT1_enriched <- WT1_enriched %>% dplyr::select(chr, start,end, name)
WT1_GR <- makeGRangesFromDataFrame(WT1_enriched, keep.extra.columns = TRUE)


#integrate BER
WT1_BER_methy_overlaps <- findOverlaps(WT1_GR, BER_hyper_merged, select = "all", ignore.strand=TRUE)
WT1_BER_methy_intersect <- GenomicRanges::intersect(WT1_GR[queryHits(WT1_BER_methy_overlaps)], BER_hyper_merged[subjectHits(WT1_BER_methy_overlaps)], ignore.strand=TRUE)

mcols(WT1_BER_methy_intersect)$name <- mcols(WT1_GR)$name[queryHits(WT1_BER_methy_overlaps)]
mcols(WT1_BER_methy_intersect)$diff.Methy <- mcols(BER_hyper_merged)$diff.Methy[queryHits(WT1_BER_methy_overlaps)]

WT1_BER_methy_intersect_df <-as.data.frame(WT1_BER_methy_intersect)
WT1_BER_methy_intersect_anno <- annotatePeak(WT1_BER_methy_intersect, TxDb = txdb, tssRegion = c(-3000,3000), annoDb = "org.Hs.eg.db")
WT1_BER_methy_intersect_anno <- as.data.frame(WT1_BER_methy_intersect_anno)

#JN
WT1_JN_methy_overlaps <- findOverlaps(WT1_GR, JN_hyper_merged, select = "all", ignore.strand=TRUE)
WT1_JN_methy_intersect <- GenomicRanges::intersect(WT1_GR[queryHits(WT1_JN_methy_overlaps)], JN_hyper_merged[subjectHits(WT1_JN_methy_overlaps)], ignore.strand=TRUE)

mcols(WT1_JN_methy_intersect)$name <- mcols(WT1_GR)$name[queryHits(WT1_JN_methy_overlaps)]
mcols(WT1_JN_methy_intersect)$diff.Methy <- mcols(JN_hyper_merged)$diff.Methy[queryHits(WT1_JN_methy_overlaps)]

WT1_JN_methy_intersect_df <-as.data.frame(WT1_JN_methy_intersect)
WT1_JN_methy_intersect_anno <- annotatePeak(WT1_JN_methy_intersect, TxDb = txdb, tssRegion = c(-3000,3000), annoDb = "org.Hs.eg.db")
WT1_JN_methy_intersect_anno <- as.data.frame(WT1_JN_methy_intersect_anno)


#SK2
WT1_SK2_methy_overlaps <- findOverlaps(WT1_GR, SK2_hyper_merged, select = "all", ignore.strand=TRUE)
WT1_SK2_methy_intersect <- GenomicRanges::intersect(WT1_GR[queryHits(WT1_SK2_methy_overlaps)], SK2_hyper_merged[subjectHits(WT1_SK2_methy_overlaps)], ignore.strand=TRUE)

mcols(WT1_SK2_methy_intersect)$name <- mcols(WT1_GR)$name[queryHits(WT1_SK2_methy_overlaps)]
mcols(WT1_SK2_methy_intersect)$diff.Methy <- mcols(SK2_hyper_merged)$diff.Methy[queryHits(WT1_SK2_methy_overlaps)]

WT1_SK2_methy_intersect_df <-as.data.frame(WT1_SK2_methy_intersect)
WT1_SK2_methy_intersect_anno <- annotatePeak(WT1_SK2_methy_intersect, TxDb = txdb, tssRegion = c(-3000,3000), annoDb = "org.Hs.eg.db")
WT1_SK2_methy_intersect_anno <- as.data.frame(WT1_SK2_methy_intersect_anno)


#tumor
WT1_Tumor291_methy_overlaps <- findOverlaps(WT1_GR, Tumor291_hyper_merged, select = "all", ignore.strand=TRUE)
WT1_Tumor291_methy_intersect <- GenomicRanges::intersect(WT1_GR[queryHits(WT1_Tumor291_methy_overlaps)], Tumor291_hyper_merged[subjectHits(WT1_Tumor291_methy_overlaps)], ignore.strand=TRUE)

mcols(WT1_Tumor291_methy_intersect)$name <- mcols(WT1_GR)$name[queryHits(WT1_Tumor291_methy_overlaps)]
mcols(WT1_Tumor291_methy_intersect)$diff.Methy <- mcols(Tumor291_hyper_merged)$diff.Methy[queryHits(WT1_Tumor291_methy_overlaps)]

WT1_Tumor291_methy_intersect_df <-as.data.frame(WT1_Tumor291_methy_intersect)
WT1_Tumor291_methy_intersect_anno <- annotatePeak(WT1_Tumor291_methy_intersect, TxDb = txdb, tssRegion = c(-3000,3000), annoDb = "org.Hs.eg.db")
WT1_Tumor291_methy_intersect_anno <- as.data.frame(WT1_Tumor291_methy_intersect_anno)

#merge
WT1_JNBER_methy <- merge(WT1_BER_methy_intersect_anno, WT1_JN_methy_intersect_anno, by="SYMBOL")
WT1_JNBERSK2_methy <- merge(WT1_JNBER_methy, WT1_SK2_methy_intersect_anno, by="SYMBOL")
WT1_JNBERSK2Tumor_methy <- merge(WT1_JNBERSK2_methy, WT1_Tumor291_methy_intersect_anno, by="SYMBOL")



#distinct genes
WT1_JNBERSK2_methy_distinct <- WT1_JNBERSK2_methy %>% distinct(SYMBOL, .keep_all = TRUE)




setwd("D:/2_RNA Seq_Files/Excel Sheets")
fold_change <- read.xlsx("Regulation_JN_BER_SK2_BOD_KTS_ALL_GENES.xlsx")
JN_fold_change <- fold_change[,c(1,2,3)]
JN_fold_change <- JN_fold_change %>% dplyr::rename(Gene = gene_name, log2_fold_change =JNshWT1_log2FC) %>%
  na.omit()


WT1_JN_methy_RNA <- merge(JN_fold_change, WT1_JN_methy_intersect_anno, by.x="Gene", by.y="SYMBOL")

#JN
#need to annotate for the gene 
setwd("D:/1_ATAC Files/ATAC 2022 Files")
JN_DEpeaksAnno <- annotatePeak("JN_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0

JN_DEpeaksAnno<-as.data.frame(JN_DEpeaksAnno) #as data frame to manipulate easier
JN_DEpeaksAnno <- JN_DEpeaksAnno %>%  #rename columns! 
  rename("FC" = "V5") #I SWEAR TO GOD THIS IS DIFFERENT EVERY TIME ON HOW TO CHANGE THE NAME

##JN_DEpeaksAnno files is the 14072 JN DAR Peaks (not genes!). want to filter by FC HERE before merging any gene stuff
JN_DEpeaks_DAR <- JN_DEpeaksAnno %>% filter(FC > 1 | FC < -1)  #filter by abs value; works great, nothing between -.99 to .99!  


WT1_JN_methy_RNA_DAR <- merge(WT1_JN_methy_RNA, JN_DEpeaks_DAR, by.x="Gene", by.y="SYMBOL")
WT1_JN_methy_RNA_DAR_distinct <- WT1_JN_methy_RNA_DAR %>% distinct(Gene, .keep_all = TRUE)
WT1_JN_methy_RNA_DAR_distinct_sig <- WT1_JN_methy_RNA_DAR_distinct %>% filter(log2_fold_change < -.99 & FC < -.99)
