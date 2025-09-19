library(openxlsx)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(VennDiagram)
library(RColorBrewer)
library(org.Hs.eg.db)
library(clusterProfiler)



##so for my methy stuff I need to make: a gene file of all names that will then be split into up and downregulated
## can probably try this in a few ways. Can try NOT filter by padj from RNA since all of my DMLs are significant. Can also try WITH filter for padj based on RNA data
##justin filtered the expression to be greater or less than |0.58|. Truly don't know where he's pulling that from or why. I think I'll just start with |1| and maybe bump down to |0.5| if needed
##I'll keep it filtered to promoter regions for now


#### Methyl load in and clean up
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
BER_hyper_merged <- bind_rows(BER_minusDox_hyper %>% dplyr::select(SYMBOL, start, end, diff.Methy, annotation, GENENAME, length, nCG),
                              BER_plusDox_hyper %>% dplyr::select(SYMBOL, start, end, diff.Methy, annotation, GENENAME, length, nCG))

JN_hyper_merged <- bind_rows(JN_minusDox_hyper %>% dplyr::select(SYMBOL, start, end, diff.Methy, annotation, GENENAME, length, nCG),
                             JN_plusDox_hyper %>% dplyr::select(SYMBOL, start, end, diff.Methy, annotation, GENENAME, length, nCG))

SK2_hyper_merged <- bind_rows(SK2_minusDox_hyper %>% dplyr::select(SYMBOL, start, end, diff.Methy, annotation, GENENAME, length, nCG),
                              SK2_plusDox_hyper %>% dplyr::select(SYMBOL, start, end, diff.Methy, annotation, GENENAME, length, nCG))

Tumor291_hyper_merged <- bind_rows(Tumor2911_hyper %>% dplyr::select(SYMBOL, start, end, diff.Methy, annotation, GENENAME, length, nCG),
                                   Tumor2912_hyper %>% dplyr::select(SYMBOL, start, end, diff.Methy, annotation, GENENAME, length, nCG))


BER_hyper_merged <- BER_hyper_merged %>%
  mutate(CpG_ID = str_c(start, end, sep = "-")) %>%
  relocate(CpG_ID, .after = SYMBOL) %>% # Relocate 'coordinates' to the 2nd position.
  dplyr::select(-start, -end) %>% # Remove the 'start' and 'end' columns.
  dplyr::rename(Gene = SYMBOL) %>% 
  dplyr::filter(grepl("Promoter", annotation)) #works like this

JN_hyper_merged <- JN_hyper_merged %>%
  mutate(CpG_ID = str_c(start, end, sep = "-")) %>%
  relocate(CpG_ID, .after = SYMBOL) %>% # Relocate 'coordinates' to the 2nd position.
  dplyr::select(-start, -end) %>% # Remove the 'start' and 'end' columns.
  dplyr::rename(Gene = SYMBOL) %>%
  dplyr::filter(grepl("Promoter", annotation)) #only in promoter regions

SK2_hyper_merged <- SK2_hyper_merged %>%
  mutate(CpG_ID = str_c(start, end, sep = "-")) %>%
  relocate(CpG_ID, .after = SYMBOL) %>% # Relocate 'coordinates' to the 2nd position.
  dplyr::select(-start, -end) %>% # Remove the 'start' and 'end' columns.
  dplyr::rename(Gene = SYMBOL) %>%
  dplyr::filter(grepl("Promoter", annotation))

Tumor291_hyper_merged <- Tumor291_hyper_merged %>%
  mutate(CpG_ID = str_c(start, end, sep = "-")) %>%
  relocate(CpG_ID, .after = SYMBOL) %>% # Relocate 'coordinates' to the 2nd position.
  dplyr::select(-start, -end) %>% # Remove the 'start' and 'end' columns.
  dplyr::rename(Gene = SYMBOL) %>%
  dplyr::filter(grepl("Promoter", annotation))


##reversing the orientation!
BER_merged_swap <- BER_hyper_merged %>%
  mutate(diff.Methy = diff.Methy * -1)

JN_merged_swap <- JN_hyper_merged %>%
  mutate(diff.Methy = diff.Methy * -1)

SK2_merged_swap <- SK2_hyper_merged %>%
  mutate(diff.Methy = diff.Methy * -1)


#load in gene name, fold change, and padj
##for BER
setwd("D:/2_RNA Seq_Files/Excel Sheets")
fold_change <- read.xlsx("Regulation_JN_BER_SK2_BOD_KTS_ALL_GENES.xlsx")

BER_fold_change <- fold_change[,c(1,6,7)]
BER_fold_change <- BER_fold_change %>% dplyr::rename(Gene = gene_name, log2_fold_change =BERshWT1_log2FC) %>%
  na.omit()

JN_fold_change <- fold_change[,c(1,3,4)]
JN_fold_change <- JN_fold_change %>% dplyr::rename(Gene = gene_name, log2_fold_change =JNshWT1_log2FC) %>%
  na.omit()


SK2_fold_change <- fold_change[,c(1,9,10)]
SK2_fold_change <- SK2_fold_change %>% dplyr::rename(Gene = gene_name, log2_fold_change =SK2shWT1_log2FC) %>%
  na.omit()


##create the different individual groups
BER_gene_up <- BER_fold_change %>% filter(log2_fold_change > 1 & BER_shWT1_padj < 0.05)
BER_gene_down <- BER_fold_change %>% filter(log2_fold_change < -1 & BER_shWT1_padj < 0.05)

BER_methy_up <- BER_merged_swap %>% filter(diff.Methy > 0.05)
BER_methy_down <- BER_merged_swap %>% filter(diff.Methy < -0.05)

BER_methy_down <- na.omit(BER_methy_down)


JN_gene_up <- JN_fold_change %>% filter(log2_fold_change > 1 & JN_shWT1_padj < 0.05)
JN_gene_down <- JN_fold_change %>% filter(log2_fold_change < -1 & JN_shWT1_padj < 0.05)

JN_methy_up <- JN_merged_swap %>% filter(diff.Methy > 0.05)
JN_methy_down <- JN_merged_swap %>% filter(diff.Methy < -0.05)

JN_methy_down <- na.omit(JN_methy_down)

SK2_gene_up <- SK2_fold_change %>% filter(log2_fold_change > 1 & SK2_shWT1_padj < 0.05)
SK2_gene_down <- SK2_fold_change %>% filter(log2_fold_change < -1 & SK2_shWT1_padj < 0.05)

SK2_methy_up <- SK2_merged_swap %>% filter(diff.Methy > 0.05)
SK2_methy_down <- SK2_merged_swap %>% filter(diff.Methy < -0.05)

SK2_methy_down <- na.omit(SK2_methy_down)


#BER
BER_RNA_methy <- merge(BER_merged_swap, BER_fold_change, by="Gene")

#BER_RNA_methy <- BER_RNA_methy %>% filter(log2_fold_change > 1 | log2_fold_change < -1) #may not need to filter before using all gene list
#BER_RNA_methy <- BER_RNA_methy %>% filter(diff.Methy > 0.05 | diff.Methy < -0.05)


##create appropriate cross of Gene UP/Methy Down to test with pathway analysis
BER_gene_up_methy_down <- merge(BER_gene_up, BER_methy_down, by="Gene")

BER_gene_up_methy_up <- merge(BER_gene_up, BER_methy_up, by="Gene")

BER_unique_gene_up_methy_down <- setdiff(BER_gene_up_methy_down$Gene, BER_gene_up_methy_up$Gene)

# Create a new dataframe with only the unique genes without overlap in Methy Up or whatever
BER_gene_up_methy_down_only <- BER_gene_up_methy_down[BER_gene_up_methy_down$Gene %in% BER_unique_gene_up_methy_down, ]

BER_gene_up_methy_down_only <- BER_gene_up_methy_down_only %>% #to check that its matching the Venn and for pathway analysis (need unique gene list)
  distinct(Gene, .keep_all = TRUE)


BER_gene_up_methy_down_genes <- as.character(BER_gene_up_methy_down_only$Gene) #pulling names/id's of filtered downregulated genes

#gene.df <- bitr(gene, fromType = "ENTREZID",
               # toType = c("ENSEMBL", "SYMBOL"),
              #  OrgDb = org.Hs.eg.db)


BER_RNA_methy_genes <- as.character(BER_RNA_methy$Gene) #all genes, from before any filtering

BER_gene_up_methy_down_ego <- enrichGO(gene = BER_gene_up_methy_down_genes, #downregulated genes (i know name is weird but just notation)
                                       universe = BER_RNA_methy_genes, #all genes, similar notation stuff
                                       keyType = "SYMBOL",
                                       OrgDb = org.Hs.eg.db, 
                                       ont = "BP", #try MF and BP
                                       pAdjustMethod = "BH", 
                                       qvalueCutoff = 0.05, 
                                       readable = TRUE)

BER_gene_up_methy_down_summary_BP <- data.frame(BER_gene_up_methy_down_ego) #just converting to df
#write.csv(BER_down_cluster_summary, "Downregulated_BER_GO_ORA.csv")
dotplot(BER_gene_up_methy_down_ego, showCategory=30, title="BER Gene Up Methy Down (MF)")



#JN
JN_RNA_methy <- merge(JN_merged_swap, JN_fold_change, by="Gene")

#JN_RNA_methy <- JN_RNA_methy %>% filter(log2_fold_change > 1 | log2_fold_change < -1) #may not need to filter before using all gene list
#JN_RNA_methy <- JN_RNA_methy %>% filter(diff.Methy > 0.05 | diff.Methy < -0.05)


##create appropriate cross of Gene UP/Methy Down to test with pathway analysis
JN_gene_up_methy_down <- merge(JN_gene_up, JN_methy_down, by="Gene")

JN_gene_up_methy_up <- merge(JN_gene_up, JN_methy_up, by="Gene")

JN_unique_gene_up_methy_down <- setdiff(JN_gene_up_methy_down$Gene, JN_gene_up_methy_up$Gene)

# Create a new dataframe with only the unique genes without overlap in Methy Up or whatever
JN_gene_up_methy_down_only <- JN_gene_up_methy_down[JN_gene_up_methy_down$Gene %in% JN_unique_gene_up_methy_down, ]

JN_gene_up_methy_down_only <- JN_gene_up_methy_down_only %>% #to check that its matching the Venn and for pathway analysis (need unique gene list)
  distinct(Gene, .keep_all = TRUE)


JN_gene_up_methy_down_genes <- as.character(JN_gene_up_methy_down_only$Gene) #pulling names/id's of filtered downregulated genes

#gene.df <- bitr(gene, fromType = "ENTREZID",
               # toType = c("ENSEMBL", "SYMBOL"),
               # OrgDb = org.Hs.eg.db)


JN_RNA_methy_genes <- as.character(JN_RNA_methy$Gene) #all genes, from before any filtering

JN_gene_up_methy_down_ego <- enrichGO(gene = JN_gene_up_methy_down_genes, #downregulated genes (i know name is weird but just notation)
                                       universe = JN_RNA_methy_genes, #all genes, similar notation stuff
                                       keyType = "SYMBOL",
                                       OrgDb = org.Hs.eg.db, 
                                       ont = "BP", #try MF and BP
                                       pAdjustMethod = "BH", 
                                       qvalueCutoff = 0.05, 
                                       readable = TRUE)

JN_gene_up_methy_down_summary_BP <- data.frame(JN_gene_up_methy_down_ego) #just converting to df
#write.csv(BER_down_cluster_summary, "Downregulated_BER_GO_ORA.csv")
dotplot(JN_gene_up_methy_down_ego, showCategory=30, title="JN: Gene Up Methy Down (MF)")



#SK2
##create overall background files by crossing RNA vs Methy
SK2_RNA_methy <- merge(SK2_merged_swap, SK2_fold_change, by="Gene")

#SK2_RNA_methy <- SK2_RNA_methy %>% filter(log2_fold_change > 1 | log2_fold_change < -1) #may not need to filter before using all gene list
#SK2_RNA_methy <- SK2_RNA_methy %>% filter(diff.Methy > 0.05 | diff.Methy < -0.05)


##create appropriate cross of Gene UP/Methy Down to test with pathway analysis
SK2_gene_up_methy_down <- merge(SK2_gene_up, SK2_methy_down, by="Gene")

SK2_gene_up_methy_up <- merge(SK2_gene_up, SK2_methy_up, by="Gene")

SK2_unique_gene_up_methy_down <- setdiff(SK2_gene_up_methy_down$Gene, SK2_gene_up_methy_up$Gene)

# Create a new dataframe with only the unique genes without overlap in Methy Up or whatever
SK2_gene_up_methy_down_only <- SK2_gene_up_methy_down[SK2_gene_up_methy_down$Gene %in% SK2_unique_gene_up_methy_down, ]

SK2_gene_up_methy_down_only <- SK2_gene_up_methy_down_only %>% #to check that its matching the Venn and for pathway analysis (need unique gene list)
  distinct(Gene, .keep_all = TRUE)


SK2_gene_up_methy_down_genes <- as.character(SK2_gene_up_methy_down_only$Gene) #pulling names/id's of filtered downregulated genes

#gene.df <- bitr(gene, fromType = "ENTREZID",
                #toType = c("ENSEMBL", "SYMBOL"),
                #OrgDb = org.Hs.eg.db)


SK2_RNA_methy_genes <- as.character(SK2_RNA_methy$Gene) #all genes, from before any filtering

SK2_gene_up_methy_down_ego <- enrichGO(gene = SK2_gene_up_methy_down_genes, #downregulated genes (i know name is weird but just notation)
                                  universe = SK2_RNA_methy_genes, #all genes, similar notation stuff
                                  keyType = "SYMBOL",
                                  OrgDb = org.Hs.eg.db, 
                                  ont = "BP", #try MF and BP
                                  pAdjustMethod = "BH", 
                                  qvalueCutoff = 0.05, 
                                  readable = TRUE)

SK2_gene_up_methy_down_summary_BP <- data.frame(SK2_gene_up_methy_down_ego) #just converting to df
#write.csv(BER_down_cluster_summary, "Downregulated_BER_GO_ORA.csv")
dotplot(SK2_gene_up_methy_down_ego, showCategory=30, title="SK2 Gene Up Methy Down (MF)")


###combining pathways from GO
write.xlsx(JN_gene_up_methy_down_summary_BP, "JN_methy_BP_summary.xlsx")
write.xlsx(BER_gene_up_methy_down_summary_BP, "BER_methy_BP_summary.xlsx")
write.xlsx(SK2_gene_up_methy_down_summary_BP, "SK2_methy_BP_summary.xlsx")


JN_BER_BP_summary <-merge(JN_gene_up_methy_down_summary_BP, BER_gene_up_methy_down_summary_BP, by="Description")
#write.xlsx(JN_BER_BP_summary, "JN_BER_BP_summary.xlsx")

SK2_JN_BP_summary <- merge(SK2_gene_up_methy_down_summary_BP, JN_gene_up_methy_down_summary_BP, by="Description")
SK2_JN_BER_BP_summary <- merge(SK2_JN_BP_summary, BER_gene_up_methy_down_summary_BP, by="Description")
write.xlsx(SK2_JN_BER_BP_summary, "JN_BER_SK2_BP_summary.xlsx")

SK2_JN_MF_summary <- merge(SK2_gene_up_methy_down_summary_MF, JN_gene_up_methy_down_summary_MF, by="Description")
SK2_JN_BER_MF_summary <- merge(SK2_JN_MF_summary, BER_gene_up_methy_down_summary_MF, by="Description")

my_colors <- brewer.pal(3, "Set1")

venn.diagram(
  x = list(SK2_gene_up_methy_down_summary_BP$Description, JN_gene_up_methy_down_summary_BP$Description, BER_gene_up_methy_down_summary_BP$Description),
  category.names = c("SK2-DSRCT", "JN-DSRCT", "BER-DSRCT"),
  filename = "JNBERSK2_ORA_GO_overlap.png",
  output = FALSE ,
  imagetype="png",
  main = "JN-BER-SK2-DSRCT ORA Overlap",
  main.cex = 1.5,
  cex = 1.5,
  cat.cex = 1.5,
  fill = my_colors,
  cat.pos = c(-20, 20, 0),
  cat.dist = c(0.05, 0.05, -0.43)
)


####KEGG 
#I ended up going back up in code above and filtering for padj < 0.05 for this stuff; can also try in Enrich GO with that too. 
##i did go back and try it 

#Prepare Enteriz gene ids
setwd("D:/2_RNA Seq_Files/Excel Sheets/DEG")
entrez_ids <- read.csv("annotations_ahb.csv")


BER_gene_up_methy_down_entrez <- merge(BER_gene_up_methy_down_only, entrez_ids, by.x="Gene", by.y="gene_name")
BER_gene_up_methy_down_entrez <- as.character(BER_gene_up_methy_down_entrez$entrezid)


## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
#can use the same one as above from EnrichGO
BER_RNA_methy_entrez <- merge(BER_RNA_methy, entrez_ids, by.x="Gene", by.y="gene_name")
BER_RNA_methy_entrez <- as.character(BER_RNA_methy_entrez$entrezid)




BER_gene_up_methy_down_KEGG <- enrichKEGG(
  gene = BER_gene_up_methy_down_entrez,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = BER_RNA_methy_entrez,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

## Output results from KEGG analysis to a table
BER_gene_up_methy_down_KEGG_summary <- data.frame(BER_gene_up_methy_down_KEGG)
#write.csv(BER_down_KEGG_cluster_summary, "BER_downregulated_genes_KEGG_analysis.csv")

## Dotplot 
dotplot(BER_gene_up_methy_down_KEGG, showCategory=20)


#JN
#Prepare Enteriz gene ids

JN_gene_up_methy_down_entrez <- merge(JN_gene_up_methy_down_only, entrez_ids, by.x="Gene", by.y="gene_name")
JN_gene_up_methy_down_entrez <- as.character(JN_gene_up_methy_down_entrez$entrezid)


## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
#can use the same one as above from EnrichGO
JN_RNA_methy_entrez <- merge(JN_RNA_methy, entrez_ids, by.x="Gene", by.y="gene_name")
JN_RNA_methy_entrez <- as.character(JN_RNA_methy_entrez$entrezid)




JN_gene_up_methy_down_KEGG <- enrichKEGG(
  gene = JN_gene_up_methy_down_entrez,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = JN_RNA_methy_entrez,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

## Output results from KEGG analysis to a table
JN_gene_up_methy_down_KEGG_summary <- data.frame(JN_gene_up_methy_down_KEGG)
#write.csv(BER_down_KEGG_cluster_summary, "BER_downregulated_genes_KEGG_analysis.csv")

## Dotplot 
dotplot(JN_gene_up_methy_down_KEGG, showCategory=20)


#SK2

SK2_gene_up_methy_down_entrez <- merge(SK2_gene_up_methy_down_only, entrez_ids, by.x="Gene", by.y="gene_name")
SK2_gene_up_methy_down_entrez <- as.character(SK2_gene_up_methy_down_entrez$entrezid)


## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
#can use the same one as above from EnrichGO
SK2_RNA_methy_entrez <- merge(SK2_RNA_methy, entrez_ids, by.x="Gene", by.y="gene_name")
SK2_RNA_methy_entrez <- as.character(SK2_RNA_methy_entrez$entrezid)




SK2_gene_up_methy_down_KEGG <- enrichKEGG(
  gene = SK2_gene_up_methy_down_entrez,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = SK2_RNA_methy_entrez,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

## Output results from KEGG analysis to a table
SK2_gene_up_methy_down_KEGG_summary <- data.frame(SK2_gene_up_methy_down_KEGG)
#write.csv(BER_down_KEGG_cluster_summary, "BER_downregulated_genes_KEGG_analysis.csv")

## Dotplot 
dotplot(SK2_gene_up_methy_down_KEGG, showCategory=20)





