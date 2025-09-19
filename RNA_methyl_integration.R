library(openxlsx)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(VennDiagram)
library(RColorBrewer)
library(org.Hs.eg.db)
library(clusterProfiler)

##loading in
setwd("D:/2_RNA Seq_Files/RSEM_TPM")
tic()
BER_ctrl1 <- read.delim("BER_shWT1_31_NT.genes.results", header = TRUE)
BER_ctrl2 <- read.delim("BER_shWT1_33_NT.genes.results", header = TRUE)
BER_test1 <- read.delim("BER_shWT1_31_4Days.genes.results", header = TRUE)
BER_test2 <- read.delim("BER_shWT1_33_4Days.genes.results", header = TRUE)

JN_ctrl1 <- read.delim("JNshWT1neg1.genes.results", header = TRUE)
JN_ctrl2 <- read.delim("JNshWT1neg2.genes.results", header = TRUE)
JN_test1 <- read.delim("JNshWT1pos1.genes.results", header = TRUE)
JN_test2 <- read.delim("JNshWT1pos2.genes.results", header = TRUE)

BOD_ctrl1 <- read.delim("BODshWT1neg1.genes.results", header = TRUE)
BOD_ctrl2 <- read.delim("BODshWT1neg2.genes.results", header = TRUE)
BOD_test1 <- read.delim("BODshWT1pos1.genes.results", header = TRUE)
BOD_test2 <- read.delim("BODshWT1pos2.genes.results", header = TRUE)

SK2_ctrl1 <- read.delim("SK2shWT1neg1.genes.results", header = TRUE)
SK2_ctrl2 <- read.delim("SK2shWT1neg2.genes.results", header = TRUE)
SK2_test1 <- read.delim("SK2shWT1pos1.genes.results", header = TRUE)
SK2_test2 <- read.delim("SK2shWT1pos2.genes.results", header = TRUE)


### removing the prefix on genes and selecting Gene and TPM values
BER_ctrl1 <- BER_ctrl1 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
BER_ctrl2 <- BER_ctrl2 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
BER_test1 <- BER_test1 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
BER_test2 <- BER_test2 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)


JN_ctrl1 <- JN_ctrl1 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
JN_ctrl2 <- JN_ctrl2 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
JN_test1 <- JN_test1 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
JN_test2 <- JN_test2 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)

BOD_ctrl1 <- BOD_ctrl1 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
BOD_ctrl2 <- BOD_ctrl2 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
BOD_test1 <- BOD_test1 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
BOD_test2 <- BOD_test2 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)


SK2_ctrl1 <- SK2_ctrl1 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
SK2_ctrl2 <- SK2_ctrl2 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
SK2_test1 <- SK2_test1 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)
SK2_test2 <- SK2_test2 %>%
  mutate(gene_id = str_remove(gene_id, "^[^_]*_")) %>%
  dplyr::select(gene_id, TPM)


##average and combine 
BER_ctrl1 <- BER_ctrl1 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

BER_ctrl2 <- BER_ctrl2 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

BER_test1 <- BER_test1 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

BER_test2 <- BER_test2 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

BER_ctrl <- full_join(BER_ctrl1, BER_ctrl2, by = "gene_id")
BER_ctrl <- BER_ctrl %>%
  mutate(Avg_TPM = rowMeans(cbind(TPM.x, TPM.y), na.rm = TRUE)) %>%
  dplyr::select(-TPM.x, -TPM.y) 

BER_test <- full_join(BER_test1, BER_test2, by = "gene_id")
BER_test <- BER_test %>%
  mutate(Avg_TPM = rowMeans(cbind(TPM.x, TPM.y), na.rm = TRUE)) %>%
  dplyr::select(-TPM.x, -TPM.y) 


JN_ctrl1 <- JN_ctrl1 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

JN_ctrl2 <- JN_ctrl2 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

JN_test1 <- JN_test1 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

JN_test2 <- JN_test2 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

JN_ctrl <- full_join(JN_ctrl1, JN_ctrl2, by = "gene_id")
JN_ctrl <- JN_ctrl %>%
  mutate(Avg_TPM = rowMeans(cbind(TPM.x, TPM.y), na.rm = TRUE)) %>%
  dplyr::select(-TPM.x, -TPM.y) 

JN_test <- full_join(JN_test1, JN_test2, by = "gene_id")
JN_test <- JN_test %>%
  mutate(Avg_TPM = rowMeans(cbind(TPM.x, TPM.y), na.rm = TRUE)) %>%
  dplyr::select(-TPM.x, -TPM.y) 


BOD_ctrl1 <- BOD_ctrl1 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

BOD_ctrl2 <- BOD_ctrl2 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

BOD_test1 <- BOD_test1 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

BOD_test2 <- BOD_test2 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

BOD_ctrl <- full_join(BOD_ctrl1, BOD_ctrl2, by = "gene_id")
BOD_ctrl <- BOD_ctrl %>%
  mutate(Avg_TPM = rowMeans(cbind(TPM.x, TPM.y), na.rm = TRUE)) %>%
  dplyr::select(-TPM.x, -TPM.y) 

BOD_test <- full_join(BOD_test1, BOD_test2, by = "gene_id")
BOD_test <- BOD_test %>%
  mutate(Avg_TPM = rowMeans(cbind(TPM.x, TPM.y), na.rm = TRUE)) %>%
  dplyr::select(-TPM.x, -TPM.y) 


SK2_ctrl1 <- SK2_ctrl1 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

SK2_ctrl2 <- SK2_ctrl2 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

SK2_test1 <- SK2_test1 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

SK2_test2 <- SK2_test2 %>%
  group_by(gene_id) %>%
  summarise(TPM = mean(TPM))

SK2_ctrl <- full_join(SK2_ctrl1, SK2_ctrl2, by = "gene_id")
SK2_ctrl <- SK2_ctrl %>%
  mutate(Avg_TPM = rowMeans(cbind(TPM.x, TPM.y), na.rm = TRUE)) %>%
  dplyr::select(-TPM.x, -TPM.y) 

SK2_test <- full_join(SK2_test1, SK2_test2, by = "gene_id")
SK2_test <- SK2_test %>%
  mutate(Avg_TPM = rowMeans(cbind(TPM.x, TPM.y), na.rm = TRUE)) %>%
  dplyr::select(-TPM.x, -TPM.y)


BER_tpm <- full_join(BER_ctrl, BER_test, by = "gene_id", suffix = c("1", "2"))
BER_tpm <- BER_tpm %>% rename("gene_id" = "Gene",
                              "Avg_TPM1" = "TPM_no_dox",
                              "Avg_TPM2" = "TPM_dox")

JN_tpm <- full_join(JN_ctrl, JN_test, by = "gene_id", suffix = c("1", "2"))
JN_tpm <- JN_tpm %>% rename("gene_id" = "Gene",
                            "Avg_TPM1" = "TPM_no_dox",
                            "Avg_TPM2" = "TPM_dox")

BOD_tpm <- full_join(BOD_ctrl, BOD_test, by = "gene_id", suffix = c("1", "2"))
BOD_tpm <- BOD_tpm %>% rename("gene_id" = "Gene",
                              "Avg_TPM1" = "TPM_no_dox",
                              "Avg_TPM2" = "TPM_dox")

SK2_tpm <- full_join(SK2_ctrl, SK2_test, by = "gene_id", suffix = c("1", "2"))
SK2_tpm <- SK2_tpm %>% rename("gene_id" = "Gene",
                              "Avg_TPM1" = "TPM_no_dox",
                              "Avg_TPM2" = "TPM_dox")

#### Methyl clean up
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





##small reset

combined_methyl_data <- rbind(
  BER_minusDox_hyper[, c("SYMBOL", "start", "end", "mean.Methy.in.Dox", "mean.Methy.in.pDox")], 
  BER_plusDox_hyper[, c("SYMBOL", "start", "end", "mean.Methy.in.Dox", "mean.Methy.in.pDox")]
)

combined_methyl_data <- combined_methyl_data %>% mutate(CpG_ID=str_c(start, end, sep = "-")) %>%
  relocate(CpG_ID, .after = SYMBOL) %>% # Relocate 'coordinates' to the 2nd position.
  dplyr::select(-start, -end) %>% # Remove the 'start' and 'end' columns.)
  dplyr::rename(Gene = SYMBOL, Methy_no_dox = mean.Methy.in.Dox, Methy_dox = mean.Methy.in.pDox) 


# 1. Merge Data (no changes needed here)
combined_data <- combined_methyl_data %>%
  left_join(BER_tpm, by = "Gene")

combined_data_complete <- combined_data %>%
  na.omit()  # Removes all rows with any NA values

combined_data_complete <- combined_data_complete[,c(2,1,3,4,5,6)]



############# making a scatter plot of diff.methy and fold change

##for BER
setwd("D:/2_RNA Seq_Files/Excel Sheets")
fold_change <- read.xlsx("Regulation_JN_BER_SK2_BOD_KTS_ALL_GENES.xlsx")
BER_fold_change <- fold_change[,c(1,6)]
BER_fold_change <- BER_fold_change %>% dplyr::rename(Gene = gene_name, log2_fold_change =BERshWT1_log2FC) %>%
  na.omit()
  

# Merge directly (no aggregation)
BER_merged_df <- BER_hyper_merged %>% 
  inner_join(BER_fold_change, by = "Gene")


# Create scatter plot
# Define the quadrant boundaries
log2_threshold <- 1
methy_threshold <- 0.1

##reversing the orientation!
BER_merged_df_swap <- BER_merged_df %>%
  mutate(diff.Methy = diff.Methy * -1)

# Create a new column for quadrant assignment
BER_merged_df_swap <- BER_merged_df_swap %>%
  mutate(Quadrant = case_when(
    log2_fold_change < -log2_threshold & diff.Methy < -methy_threshold ~ "Q1",  # Bottom left
    log2_fold_change > log2_threshold  & diff.Methy < -methy_threshold ~ "Q2",  # Top left
    log2_fold_change < -log2_threshold & diff.Methy > methy_threshold  ~ "Q3",  # bottom right
    log2_fold_change > log2_threshold  & diff.Methy > methy_threshold  ~ "Q4",  # Top right
    TRUE                                                             ~ "NS"   # Not significant
  ))

BER_merged_df_swap <- BER_merged_df_swap %>%
  mutate(Methylation = case_when(diff.Methy > 0 ~ "Hyper",
                                 diff.Methy <0 ~ "Hypo"))
# Find top 5 points within each quadrant
BER_top_points <- BER_merged_df_swap %>%
  group_by(Quadrant) %>%
  arrange(desc(abs(diff.Methy))) %>%  # Sort by absolute Diff.methy, descending
  slice_head(n = 10)  # Take top 5 points in each group


ggplot(BER_merged_df_swap, aes(x = diff.Methy, y = log2_fold_change, color = Quadrant)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Q1" = "blue", "Q2" = "red", "Q3" = "#006400", "Q4" = "purple", "NS" = "gray"),
                     labels = c("Q1" = "Hypo-Down", 
                                "Q2" = "Hypo-Up", 
                                "Q3" = "Hyper-Down",
                                "Q4" = "Hyper-Up",
                                "NS" = "Not Significant")) +  # Custom colors
  labs(title = "BER: Methylation Difference vs. Log2 Expression Change after EWSR1-WT1 KD",
       x = "Methylation Difference",
       y = "Log2 Expression Change") +
  theme_minimal() +
  geom_hline(yintercept = c(-log2_threshold, log2_threshold), linetype = "dashed", color = "gray") +  # Add horizontal lines
  geom_vline(xintercept = c(-methy_threshold, methy_threshold), linetype = "dashed", color = "gray") +  # Add vertical lines
  geom_text_repel(data = BER_top_points, aes(label = Gene), size = 3) +
  annotate("text", x = -0.45, y = 5., label = "Hypo-Up", color = "black") +  # Adjust coordinates as needed
  annotate("text", x = 0.45, y = 5.3, label = "Hyper-Up", color = "black") +    # Adjust coordinates as needed
  annotate("text", x = -0.45, y = -5.3, label = "Hypo-Down", color = "black") +  # Adjust coordinates as needed
  annotate("text", x = 0.45, y = -5.3, label = "Hyper-Down", color = "black") +    # Adjust coordinates as needed
  theme(legend.position = "none")  # Remove the legend

BER_hyper_down <- BER_merged_df_swap %>% filter(Quadrant == "Q3")
BER_hyper_down_consensus_downreg <- merge(BER_hyper_down, consensus_downregulated, by="Gene")

BER_hypo_up <- BER_merged_df_swap %>% filter(Quadrant =="Q2")
BER_hypo_up_consensus_upreg <- merge(BER_hypo_up, consensus_upregulated, by="Gene")

###for JN
setwd("D:/2_RNA Seq_Files/Excel Sheets")
fold_change <- read.xlsx("Regulation_JN_BER_SK2_BOD_KTS_ALL_GENES.xlsx")
JN_fold_change <- fold_change[,c(1,3)]
JN_fold_change <- JN_fold_change %>% dplyr::rename(Gene = gene_name, log2_fold_change =JNshWT1_log2FC) %>%
  na.omit()


# Merge directly (no aggregation)
JN_merged_df <- JN_hyper_merged %>% 
  inner_join(JN_fold_change, by = "Gene")


#high_DMR <- merged_df %>% filter(abs(log2_fold_change) > 1, abs(diff.Methy) > 0.2)
#higher_DMR <- merged_df %>% filter(abs(log2_fold_change) >1.5, abs(diff.Methy) >0.3)

# Create scatter plot
# Define the quadrant boundaries
log2_threshold <- 1
methy_threshold <- 0.1

##reversing the orientation!
JN_merged_df_swap <- JN_merged_df %>%
  mutate(diff.Methy = diff.Methy * -1)

# Create a new column for quadrant assignment
JN_merged_df_swap <- JN_merged_df_swap %>%
  mutate(Quadrant = case_when(
    log2_fold_change < -log2_threshold & diff.Methy < -methy_threshold ~ "Q1",  # Bottom left
    log2_fold_change > log2_threshold  & diff.Methy < -methy_threshold ~ "Q2",  # top left
    log2_fold_change < -log2_threshold & diff.Methy > methy_threshold  ~ "Q3",  # bottom right
    log2_fold_change > log2_threshold  & diff.Methy > methy_threshold  ~ "Q4",  # Top right
    TRUE                                                             ~ "NS"   # Not significant
  ))


# Find top points within each quadrant
JN_top_points <- JN_merged_df_swap %>%
  group_by(Quadrant) %>%
  arrange(desc(abs(diff.Methy))) %>%  # Sort by absolute Diff.methy, descending
  slice_head(n = 10)  # Take top 5 points in each group


ggplot(JN_merged_df_swap, aes(x = diff.Methy, y = log2_fold_change, color = Quadrant)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Q1" = "blue", "Q2" = "red", "Q3" = "#006400", "Q4" = "purple", "NS" = "gray"),
                     labels = c("Q1" = "Hypo-Down", 
                                "Q2" = "Hypo-Up", 
                                "Q3" = "Hyper-Down",
                                "Q4" = "Hyper-Up",
                                "NS" = "Not Significant")) +  # Custom colors
  labs(title = "JN: Methylation Difference vs. Log2 Expression Change after EWSR1-WT1 KD",
       x = "Methylation Difference",
       y = "Log2 Expression Change") +
  theme_minimal() +
  geom_hline(yintercept = c(-log2_threshold, log2_threshold), linetype = "dashed", color = "gray") +  # Add horizontal lines
  geom_vline(xintercept = c(-methy_threshold, methy_threshold), linetype = "dashed", color = "gray") +  # Add vertical lines
  geom_text_repel(data = JN_top_points, aes(label = Gene), size = 3) +
  annotate("text", x = -0.5, y = 5.5, label = "Hypo-Up", color = "black") +  # Adjust coordinates as needed
  annotate("text", x = 0.5, y = 5.5, label = "Hyper-Up", color = "black") +    # Adjust coordinates as needed
  annotate("text", x = -0.5, y = -5.5, label = "Hypo-Down", color = "black") +  # Adjust coordinates as needed
  annotate("text", x = 0.5, y = -5.5, label = "Hyper-Down", color = "black") +    # Adjust coordinates as needed
  theme(legend.position = "none")  # Remove the legend

JN_hypo_up <- JN_merged_df_swap %>% filter(Quadrant == "Q2")
JN_hypo_up_consensus_upreg <- merge(JN_hypo_up, consensus_upregulated, by="Gene")


##for SK2
setwd("D:/2_RNA Seq_Files/Excel Sheets")
fold_change <- read.xlsx("Regulation_JN_BER_SK2_BOD_KTS_ALL_GENES.xlsx")
SK2_fold_change <- fold_change[,c(1,9)]
SK2_fold_change <- SK2_fold_change %>% dplyr::rename(Gene = gene_name, log2_fold_change =SK2shWT1_log2FC) %>%
  na.omit()


# Merge directly (no aggregation)
SK2_merged_df <- SK2_hyper_merged %>% 
  inner_join(SK2_fold_change, by = "Gene")


# Create scatter plot
# Define the quadrant boundaries
log2_threshold <- 1
methy_threshold <- 0.1

##reversing the orientation!
SK2_merged_df_swap <- SK2_merged_df %>%
  mutate(diff.Methy = diff.Methy * -1)


# Create a new column for quadrant assignment
SK2_merged_df_swap <- SK2_merged_df_swap %>%
  mutate(Quadrant = case_when(
    log2_fold_change < -log2_threshold & diff.Methy < -methy_threshold ~ "Q1",  # Bottom left
    log2_fold_change > log2_threshold  & diff.Methy < -methy_threshold ~ "Q2",  # top left
    log2_fold_change < -log2_threshold & diff.Methy > methy_threshold  ~ "Q3",  # bottom right
    log2_fold_change > log2_threshold  & diff.Methy > methy_threshold  ~ "Q4",  # Top right
    TRUE                                                             ~ "NS"   # Not significant
  ))


# Find top 5 points within each quadrant
SK2_top_points <- SK2_merged_df_swap %>%
  group_by(Quadrant) %>%
  arrange(desc(abs(diff.Methy))) %>%  # Sort by absolute Diff.methy, descending
  slice_head(n = 10)  # Take top 5 points in each group


ggplot(SK2_merged_df_swap, aes(x = diff.Methy, y = log2_fold_change, color = Quadrant)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Q1" = "blue", "Q2" = "red", "Q3" = "#006400", "Q4" = "purple", "NS" = "gray"),
                     labels = c("Q1" = "Hypo-Down", 
                                "Q2" = "Hypo-Up", 
                                "Q3" = "Hyper-Down",
                                "Q4" = "Hyper-Up",
                                "NS" = "Not Significant")) +  # Custom colors
  labs(title = "SK2: Methylation Difference vs. Log2 Expression Change after EWSR1-WT1 KD",
       x = "Methylation Difference",
       y = "Log2 Expression Change") +
  theme_minimal() +
  geom_hline(yintercept = c(-log2_threshold, log2_threshold), linetype = "dashed", color = "gray") +  # Add horizontal lines
  geom_vline(xintercept = c(-methy_threshold, methy_threshold), linetype = "dashed", color = "gray") +  # Add vertical lines
  geom_text_repel(data = SK2_top_points, aes(label = Gene), size = 3) +
  annotate("text", x = -0.5, y = 5.5, label = "Hypo-Up", color = "black") +  # Adjust coordinates as needed
  annotate("text", x = 0.5, y = 5.5, label = "Hyper-Up", color = "black") +    # Adjust coordinates as needed
  annotate("text", x = -0.5, y = -5.5, label = "Hypo-Down", color = "black") +  # Adjust coordinates as needed
  annotate("text", x = 0.5, y = -5.5, label = "Hyper-Down", color = "black") +    # Adjust coordinates as needed
  theme(legend.position = "none")  # Remove the legend



SK2_hyper_down <- SK2_merged_df_swap %>% filter(Quadrant == "Q3")
SK2_hyper_down_consensus_downreg <- merge(SK2_hyper_down, consensus_downregulated, by="Gene")

SK2_hypo_up <- SK2_merged_df_swap %>% filter(Quadrant == "Q2")
SK2_hypo_up_consensus_upreg <- merge(SK2_hypo_up, consensus_upregulated, by="Gene")

###for Tumor 291
setwd("D:/2_RNA Seq_Files/Excel Sheets/DEG")
Tumor291_fold_change <- read.xlsx("291_v_LP9_rsem_DGE_results.xlsx")
Tumor291_fold_change <- Tumor291_fold_change[,c(1,3)]
Tumor291_fold_change <- Tumor291_fold_change %>% dplyr::rename(Gene = gene_name, log2_fold_change = log2FoldChange) %>%
  na.omit()


# Merge directly (no aggregation)
Tumor291_merged_df <- Tumor291_hyper_merged %>% 
  inner_join(Tumor291_fold_change, by = "Gene")


# Create scatter plot
# Define the quadrant boundaries
log2_threshold <- 1
methy_threshold <- 0.1

##reversing the orientation!
Tumor291_merged_df_swap <- Tumor291_merged_df %>%
  mutate(diff.Methy = diff.Methy * -1)


# Create a new column for quadrant assignment
Tumor291_merged_df_swap <- Tumor291_merged_df_swap %>%
  mutate(Quadrant = case_when(
    log2_fold_change < -log2_threshold & diff.Methy < -methy_threshold ~ "Q1",  # Bottom left
    log2_fold_change > log2_threshold  & diff.Methy < -methy_threshold ~ "Q2",  # top left
    log2_fold_change < -log2_threshold & diff.Methy > methy_threshold  ~ "Q3",  # bottom right
    log2_fold_change > log2_threshold  & diff.Methy > methy_threshold  ~ "Q4",  # Top right
    TRUE                                                             ~ "NS"   # Not significant
  ))


# Find top 5 points within each quadrant
Tumor291_top_points <- Tumor291_merged_df_swap %>%
  group_by(Quadrant) %>%
  arrange(desc(abs(diff.Methy))) %>%  # Sort by absolute Diff.methy, descending
  slice_head(n = 10)  # Take top 5 points in each group


ggplot(Tumor291_merged_df_swap, aes(x = diff.Methy, y = log2_fold_change, color = Quadrant)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Q1" = "blue", "Q2" = "red", "Q3" = "#006400", "Q4" = "purple", "NS" = "gray"),
                     labels = c("Q1" = "Hypo-Down", 
                                "Q2" = "Hypo-Up", 
                                "Q3" = "Hyper-Down",
                                "Q4" = "Hyper-Up",
                                "NS" = "Not Significant")) +  # Custom colors
  labs(title = "Tumor291: Methylation Difference vs. Log2 Expression Change in Primary vs Metastatic Tumors",
       x = "Methylation Difference",
       y = "Log2 Expression Change") +
  theme_minimal() +
  geom_hline(yintercept = c(-log2_threshold, log2_threshold), linetype = "dashed", color = "gray") +  # Add horizontal lines
  geom_vline(xintercept = c(-methy_threshold, methy_threshold), linetype = "dashed", color = "gray") +  # Add vertical lines
  geom_text_repel(data = Tumor291_top_points, aes(label = Gene), size = 3) +
  annotate("text", x = -0.75, y = 10, label = "Hypo-Up", color = "black") +  # Adjust coordinates as needed
  annotate("text", x = 0.75, y = 10, label = "Hyper-Up", color = "black") +    # Adjust coordinates as needed
  annotate("text", x = -0.75, y = -10, label = "Hypo-Down", color = "black") +  # Adjust coordinates as needed
  annotate("text", x = 0.75, y = -10, label = "Hyper-Down", color = "black") +    # Adjust coordinates as needed
  theme(legend.position = "none")  # Remove the legend



Tumor291_hyper_down <- Tumor291_merged_df_swap %>% filter(Quadrant == "Q3") 

Tumor291_hyper_down_consensus_downreg <- merge(Tumor291_hyper_down, consensus_downregulated, by="Gene")
Tumor291_hyper_down_consensus_upreg  <- merge(Tumor291_hyper_down, consensus_upregulated, by="Gene")




##
JNBER_hypo_up <- merge(JN_hypo_up_consensus_upreg, BER_hypo_up_consensus_upreg, by="Gene")
JNBERSK2_hypo_up_consensus_upreg <- merge(JNBER_hypo_up, SK2_hypo_up_consensus_upreg, by="Gene")

##
##dividing into quadrants and doing pathway analysis. Will pull out all the relevant info and break down from there
#Q1 = "Downregulated & Hypomethylated
#Q2 = "Upregulated & Hypomethylated", 
#Q3 = "Downregulated & Hypermethylated",
#Q4 = "Upregulated & Hypermethylated"

setwd("C:/Users/Flinchum/Documents/Tulane Docs/Lee Lab/Results/Survival Analysis")
survival <- read.xlsx("survival_analysis_all_genes.xlsx")

BER_merged_df_swap_surv <- merge(BER_merged_df_swap, survival, by="Gene")
BER_Q1 <- BER_merged_df_swap_surv %>%
  filter(Quadrant == "Q1") 

BER_Q2 <- BER_merged_df_swap_surv %>%
  filter(Quadrant == "Q2")

BER_Q3 <- BER_merged_df_swap_surv %>%
  filter(Quadrant == "Q3") 

BER_Q4 <- BER_merged_df_swap_surv %>%
  filter(Quadrant == "Q4")

JN_merged_df_swap_surv <- merge(JN_merged_df_swap, survival, by="Gene")
JN_Q1 <- JN_merged_df_swap_surv %>%
  filter(Quadrant == "Q1") 

JN_Q2 <- JN_merged_df_swap_surv %>%
  filter(Quadrant == "Q2")

JN_Q3 <- JN_merged_df_swap_surv %>%
  filter(Quadrant == "Q3") 

JN_Q4 <- JN_merged_df_swap_surv %>%
  filter(Quadrant == "Q4")

SK2_merged_df_swap_surv <- merge(SK2_merged_df_swap, survival, by="Gene")
SK2_Q1 <- SK2_merged_df_swap_surv %>%
  filter(Quadrant == "Q1") 

SK2_Q2 <- SK2_merged_df_swap_surv %>%
  filter(Quadrant == "Q2")

SK2_Q3 <- SK2_merged_df_swap_surv %>%
  filter(Quadrant == "Q3") 

SK2_Q4 <- SK2_merged_df_swap_surv %>%
  filter(Quadrant == "Q4")

Tumor291_merged_df_swap_surv <- merge(Tumor291_merged_df_swap, survival, by="Gene")
Tumor291_Q1 <- Tumor291_merged_df_swap_surv %>%
  filter(Quadrant == "Q1") 

Tumor291_Q2 <- Tumor291_merged_df_swap_surv %>%
  filter(Quadrant == "Q2")

Tumor291_Q3 <- Tumor291_merged_df_swap_surv %>%
  filter(Quadrant == "Q3") 

Tumor291_Q4 <- Tumor291_merged_df_swap_surv %>%
  filter(Quadrant == "Q4")

####
##trying GSEA with the Q3 stuff/maybe other Qs
## Merge the AnnotationHub dataframe with the results 
setwd("D:/2_RNA Seq_Files/Excel Sheets/DEG")
#Tumor_v_LP9 <- read.xlsx("291_v_LP9_rsem_DGE_results.xlsx")
annotations_ahb <-read.csv("annotations_ahb.csv")
annotations_ahb <- annotations_ahb %>% dplyr::rename("Gene"="gene_name")
#Tumor_GSEA <- merge(Tumor_v_LP9,annotations_ahb,by="gene_name")

Tumor291_Q3_gsea <- merge(Tumor291_merged_df_swap, annotations_ahb, by="Gene")


Tumor291_Q3_gsea <- Tumor291_Q3_gsea %>% filter(diff.Methy > 0)

Tumor291_Q3_gsea <- Tumor291_Q3_gsea[which(duplicated(Tumor291_Q3_gsea$entrezid) == F), ]


Tumor291_Q3_log2fc <- as.numeric(Tumor291_Q3_gsea$log2_fold_change)
names(Tumor291_Q3_log2fc) <- as.vector(Tumor291_Q3_gsea$entrezid)


Tumor291_Q3_log2fc <-sort(Tumor291_Q3_log2fc, decreasing= TRUE)
head(Tumor291_Q3_log2fc)

## GSEA using gene sets from KEGG pathways
Tumor291_gseaKEGG <- gseKEGG(geneList =Tumor291_Q3_log2fc, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    #nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)


## Extract the GSEA results
Tumor291_gseaKEGG_results <- Tumor291_gseaKEGG@result

#write.csv(Tumor_gseaKEGG_results, "291Tumor_v_LP9_gsea_kegg.csv", quote=F)

#Plot Some Results
#gseaplot(Tumor291_gseaKEGG, geneSetID = 'hsa04820')
dotplot(Tumor291_gseaKEGG, showCategory=20, title="Tumor291 KEGG")



# GSEA using gene sets associated with BP Gene Ontology terms
Tumor291_gseaGO <- gseGO(geneList = Tumor291_Q3_log2fc, 
                OrgDb = org.Hs.eg.db, 
                ont = 'MF', 
                #nPerm = 1000, 
                minGSSize = 20, 
                pvalueCutoff = 0.05,
                verbose = FALSE) 

Tumor291_gseaGO_MF_results <- Tumor291_gseaGO@result

write.csv(Tumor291_gseaGO_MF_results, "Tumor291_v_LP9_GO_GSEA.csv", quote=F)

#gseaplot(Tumor291_gseaGO, geneSetID = 'GO:0098609')
dotplot(Tumor291_gseaGO, showCategory=20, title = "Tumor291 GO MF Results")


# Split the core_enrichment string into individual Entrez IDs
Tumor291_gseaGO_MF_results$entrez_ids <- strsplit(Tumor291_gseaGO_MF_results$core_enrichment, "/")

# Convert the list of Entrez IDs into a vector
Tumor291_gseaGO_MF_results$entrez_ids <- lapply(Tumor291_gseaGO_MF_results$entrez_ids, unlist)

# Convert Entrez IDs to gene symbols
Tumor291_gseaGO_MF_results$gene_symbols <- lapply(Tumor291_gseaGO_MF_results$entrez_ids, function(ids) {
  mapIds(org.Hs.eg.db, keys=ids, column="SYMBOL", keytype="ENTREZID", multiVals="first")
})

Tumor291_GOMF_genes <- Tumor291_gseaGO_MF_results [,c(1,2,11)]
Tumor291_GOMF_genes <- Tumor291_GOMF_genes %>% rename("ID"="GO_ID", "core_enrichment"="Gene") %>%
  separate_rows(Gene, sep = "/") # Assuming '/' separates genes

Tumor291_GOMF_genes$Gene <- lapply(Tumor291_GOMF_genes$Gene, unlist)
tic("hi test")
Tumor291_GOMF_genes$Gene <- lapply(Tumor291_GOMF_genes$Gene, function(ids) {
  mapIds(org.Hs.eg.db, keys=ids, column="SYMBOL", keytype="ENTREZID", multiVals="first")
})
toc()

##pretty sure each gene symbol in the "Gene" column was its own list with one entry. Had to sort out before being able to use
# Unlist the 'Gene' column and create a new data frame
Tumor291_MF_unlisted_genes <- unlist(Tumor291_GOMF_genes$Gene)
#Tumor291_gene_symbols <- data.frame(Gene = Tumor291_unlisted_genes) #basically just to inspect/check if needed

# Create indices to repeat the GO_IDs and Descriptions
Tumor291_MF_indices <- rep(1:nrow(Tumor291_GOMF_genes), lengths(Tumor291_GOMF_genes$Gene))

# Repeat GO_IDs and Descriptions to match the unlisted genes
Tumor291_GOMF_genes <- data.frame(
  GO_ID = Tumor291_GOMF_genes$GO_ID[Tumor291_MF_indices],
  Description = Tumor291_GOMF_genes$Description[Tumor291_MF_indices],
  Gene = Tumor291_MF_unlisted_genes
)

Tumor291_GOMF_genes_surv <- merge(Tumor291_GOMF_genes, survival, by="Gene")

##BER gsea
#setwd("D:/2_RNA Seq FIles/Excel Sheets/DEG")
#annotations_ahb <-read.csv("annotations_ahb.csv")
#annotations_ahb <- annotations_ahb %>% rename("gene_name"="Gene")

BER_Q3_gsea <- merge(BER_merged_df_swap, annotations_ahb, by="Gene")

BER_Q3_gsea <- BER_Q3_gsea %>% filter(diff.Methy > 0)

#Remove Duplicated Genes Method 1 - Removes way more "duplicates" than method 2
BER_Q3_gsea <- BER_Q3_gsea[which(duplicated(BER_Q3_gsea$entrezid) == F), ]

BER_Q3_log2fc <- as.numeric(BER_Q3_gsea$log2_fold_change)
names(BER_Q3_log2fc) <- as.vector(BER_Q3_gsea$entrezid)


#Sort the list by log2foldchange
BER_Q3_log2fc <-sort(BER_Q3_log2fc, decreasing= TRUE)
head(BER_Q3_log2fc)


BER_gseaKEGG <- gseKEGG(geneList =BER_Q3_log2fc, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    #nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
BER_gseaKEGG_results <- BER_gseaKEGG@result

#write.csv(Tumor_gseaKEGG_results, "291Tumor_v_LP9_gsea_kegg.csv", quote=F)

#Plot Some Results
#gseaplot(BER_gseaKEGG, geneSetID = 'hsa05206')
dotplot(BER_gseaKEGG, showCategory=30, title="BER KEGG")

# GSEA using gene sets associated with BP Gene Ontology terms
BER_gseaGO <- gseGO(geneList = BER_Q3_log2fc, 
                         OrgDb = org.Hs.eg.db, 
                         ont = 'BP', 
                         #nPerm = 1000, 
                         minGSSize = 20, 
                         pvalueCutoff = 0.05,
                         verbose = FALSE) 

BER_gseaGO_results <- BER_gseaGO@result

write.xlsx(BER_gseaGO_results, "BER_methy_GO_GSEA.xlsx", quote=F)

gseaplot(BER_gseaGO, geneSetID = 'GO:0048762')
dotplot(BER_gseaGO, showCategory=30, title = "BER GO Results")


# Split the core_enrichment string into individual Entrez IDs
BER_gseaGO_results$entrez_ids <- strsplit(BER_gseaGO_results$core_enrichment, "/")

# Convert the list of Entrez IDs into a vector
BER_gseaGO_results$entrez_ids <- lapply(BER_gseaGO_results$entrez_ids, unlist)

# Convert Entrez IDs to gene symbols
BER_gseaGO_results$gene_symbols <- lapply(BER_gseaGO_results$entrez_ids, function(ids) {
  mapIds(org.Hs.eg.db, keys=ids, column="SYMBOL", keytype="ENTREZID", multiVals="first")
})

BER_GO_genes <- BER_gseaGO_results [,c(1,2,11)]
BER_GO_genes <- BER_GO_genes %>% rename("ID"="GO_ID", "core_enrichment"="Gene") %>%
  separate_rows(Gene, sep = "/") # Assuming '/' separates genes

BER_GO_genes$Gene <- lapply(BER_GO_genes$Gene, unlist)
BER_GO_genes$Gene <- lapply(BER_GO_genes$Gene, function(ids) {
  mapIds(org.Hs.eg.db, keys=ids, column="SYMBOL", keytype="ENTREZID", multiVals="first")
})


##pretty sure each gene symbol in the "Gene" column was its own list with one entry. Had to sort out before being able to use
# Unlist the 'Gene' column and create a new data frame
BER_unlisted_genes <- unlist(BER_GO_genes$Gene)
#BER_gene_symbols <- data.frame(Gene = BER_unlisted_genes) #basically just to inspect/check if needed

# Create indices to repeat the GO_IDs and Descriptions
BER_indices <- rep(1:nrow(BER_GO_genes), lengths(BER_GO_genes$Gene))

# Repeat GO_IDs and Descriptions to match the unlisted genes
BER_GO_genes <- data.frame(
  GO_ID = BER_GO_genes$GO_ID[BER_indices],
  Description = BER_GO_genes$Description[BER_indices],
  Gene = BER_unlisted_genes
)




##SK2 gsea
#setwd("D:/2_RNA Seq FIles/Excel Sheets/DEG")
#annotations_ahb <-read.csv("annotations_ahb.csv")
#annotations_ahb <- annotations_ahb %>% rename("gene_name"="Gene")

SK2_Q3_gsea <- merge(SK2_merged_df_swap, annotations_ahb, by="Gene")

SK2_Q3_gsea <- SK2_Q3_gsea %>% filter(diff.Methy > 0)

#Remove Duplicated Genes Method 1 - Removes way more "duplicates" than method 2
SK2_Q3_gsea <- SK2_Q3_gsea[which(duplicated(SK2_Q3_gsea$entrezid) == F), ]

SK2_Q3_log2fc <- as.numeric(SK2_Q3_gsea$log2_fold_change)
names(SK2_Q3_log2fc) <- as.vector(SK2_Q3_gsea$entrezid)


#Sort the list by log2foldchange
SK2_Q3_log2fc <-sort(SK2_Q3_log2fc, decreasing= TRUE)
head(SK2_Q3_log2fc)


SK2_gseaKEGG <- gseKEGG(geneList =SK2_Q3_log2fc, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    #nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
SK2_gseaKEGG_results <- SK2_gseaKEGG@result

#write.csv(Tumor_gseaKEGG_results, "291Tumor_v_LP9_gsea_kegg.csv", quote=F)

#Plot Some Results
#gseaplot(SK2_gseaKEGG, geneSetID = 'hsa05206')
dotplot(SK2_gseaKEGG, showCategory=30, title="SK2 KEGG")


# GSEA using gene sets associated with BP Gene Ontology terms
SK2_gseaGO <- gseGO(geneList = SK2_Q3_log2fc, 
                    OrgDb = org.Hs.eg.db, 
                    ont = 'BP', 
                    #nPerm = 1000, 
                    minGSSize = 20, 
                    pvalueCutoff = 0.05,
                    verbose = FALSE) 

SK2_gseaGO_results <- SK2_gseaGO@result

write.csv(SK2_gseaGO_results, "SK2_v_LP9_GO_GSEA.csv", quote=F)

#gseaplot(SK2_gseaGO, geneSetID = 'GO:0098609')
dotplot(SK2_gseaGO, showCategory=30, title = "SK2 GO Results")


# Split the core_enrichment string into individual Entrez IDs
SK2_gseaGO_results$entrez_ids <- strsplit(SK2_gseaGO_results$core_enrichment, "/")

# Convert the list of Entrez IDs into a vector
SK2_gseaGO_results$entrez_ids <- lapply(SK2_gseaGO_results$entrez_ids, unlist)

# Convert Entrez IDs to gene symbols
SK2_gseaGO_results$gene_symbols <- lapply(SK2_gseaGO_results$entrez_ids, function(ids) {
  mapIds(org.Hs.eg.db, keys=ids, column="SYMBOL", keytype="ENTREZID", multiVals="first")
})

SK2_GO_genes <- SK2_gseaGO_results [,c(1,2,11)]
SK2_GO_genes <- SK2_GO_genes %>% rename("ID"="GO_ID", "core_enrichment"="Gene") %>%
  separate_rows(Gene, sep = "/") # Assuming '/' separates genes

SK2_GO_genes$Gene <- lapply(SK2_GO_genes$Gene, unlist)
SK2_GO_genes$Gene <- lapply(SK2_GO_genes$Gene, function(ids) {
  mapIds(org.Hs.eg.db, keys=ids, column="SYMBOL", keytype="ENTREZID", multiVals="first")
})


##pretty sure each gene symbol in the "Gene" column was its own list with one entry. Had to sort out before being able to use
# Unlist the 'Gene' column and create a new data frame
SK2_unlisted_genes <- unlist(SK2_GO_genes$Gene)
#SK2_gene_symbols <- data.frame(Gene = SK2_unlisted_genes) #basically just to inspect/check if needed

# Create indices to repeat the GO_IDs and Descriptions
SK2_indices <- rep(1:nrow(SK2_GO_genes), lengths(SK2_GO_genes$Gene))

# Repeat GO_IDs and Descriptions to match the unlisted genes
SK2_GO_genes <- data.frame(
  GO_ID = SK2_GO_genes$GO_ID[SK2_indices],
  Description = SK2_GO_genes$Description[SK2_indices],
  Gene = SK2_unlisted_genes
)




##JN gsea
#setwd("D:/2_RNA Seq FIles/Excel Sheets/DEG")
#annotations_ahb <-read.csv("annotations_ahb.csv")
#annotations_ahb <- annotations_ahb %>% rename("gene_name"="Gene")

JN_Q3_gsea <- merge(JN_merged_df_swap, annotations_ahb, by="Gene")

JN_Q3_gsea <- JN_Q3_gsea %>% filter(diff.Methy > 0)


#Remove Duplicated Genes Method 1 - Removes way more "duplicates" than method 2
JN_Q3_gsea <- JN_Q3_gsea[which(duplicated(JN_Q3_gsea$entrezid) == F), ]

JN_Q3_log2fc <- as.numeric(JN_Q3_gsea$log2_fold_change)
names(JN_Q3_log2fc) <- as.vector(JN_Q3_gsea$entrezid)


#Sort the list by log2foldchange
JN_Q3_log2fc <-sort(JN_Q3_log2fc, decreasing= TRUE)
head(JN_Q3_log2fc)


JN_gseaKEGG <- gseKEGG(geneList =JN_Q3_log2fc, # ordered named vector of fold changes (Entrez IDs are the associated names)
                        organism = "hsa", # supported organisms listed below
                        #nPerm = 1000, # default number permutations
                        minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                        pvalueCutoff = 0.05, # padj cutoff value
                        verbose = FALSE)

## Extract the GSEA results
JN_gseaKEGG_results <- JN_gseaKEGG@result

#write.csv(Tumor_gseaKEGG_results, "291Tumor_v_LP9_gsea_kegg.csv", quote=F)

#Plot Some Results
#gseaplot(JN_gseaKEGG, geneSetID = 'hsa05206')
dotplot(JN_gseaKEGG, showCategory=30, title="JN KEGG")


# GSEA using gene sets associated with BP Gene Ontology terms
JN_gseaGO <- gseGO(geneList = JN_Q3_log2fc, 
                    OrgDb = org.Hs.eg.db, 
                    ont = 'BP', 
                    #nPerm = 1000, 
                    minGSSize = 20, 
                    pvalueCutoff = 0.05,
                    verbose = FALSE) 

JN_gseaGO_results <- JN_gseaGO@result

#write.csv(Tumor_gseaGO_results, "JN_v_LP9_GO_GSEA.csv", quote=F)

#gseaplot(JN_gseaGO, geneSetID = 'GO:0098609')
dotplot(JN_gseaGO, showCategory=30, title = "JN GO Results")


# Split the core_enrichment string into individual Entrez IDs
JN_gseaGO_results$entrez_ids <- strsplit(JN_gseaGO_results$core_enrichment, "/")

# Convert the list of Entrez IDs into a vector
JN_gseaGO_results$entrez_ids <- lapply(JN_gseaGO_results$entrez_ids, unlist)

# Convert Entrez IDs to gene symbols
JN_gseaGO_results$gene_symbols <- lapply(JN_gseaGO_results$entrez_ids, function(ids) {
  mapIds(org.Hs.eg.db, keys=ids, column="SYMBOL", keytype="ENTREZID", multiVals="first")
})

JN_GO_genes <- JN_gseaGO_results [,c(1,2,11)]
JN_GO_genes <- JN_GO_genes %>% rename("ID"="GO_ID", "core_enrichment"="Gene") %>%
  separate_rows(Gene, sep = "/") # Assuming '/' separates genes

JN_GO_genes$Gene <- lapply(JN_GO_genes$Gene, unlist)
JN_GO_genes$Gene <- lapply(JN_GO_genes$Gene, function(ids) {
  mapIds(org.Hs.eg.db, keys=ids, column="SYMBOL", keytype="ENTREZID", multiVals="first")
})


##pretty sure each gene symbol in the "Gene" column was its own list with one entry. Had to sort out before being able to use
# Unlist the 'Gene' column and create a new data frame
JN_unlisted_genes <- unlist(JN_GO_genes$Gene)
#JN_gene_symbols <- data.frame(Gene = JN_unlisted_genes) #basically just to inspect/check if needed

# Create indices to repeat the GO_IDs and Descriptions
JN_indices <- rep(1:nrow(JN_GO_genes), lengths(JN_GO_genes$Gene))

# Repeat GO_IDs and Descriptions to match the unlisted genes
JN_GO_genes <- data.frame(
  GO_ID = JN_GO_genes$GO_ID[JN_indices],
  Description = JN_GO_genes$Description[JN_indices],
  Gene = JN_unlisted_genes
)



## merging the GO stuff together 
BERSK2_GO_Genes <- merge(BER_GO_genes, SK2_GO_genes, by="Gene", all=TRUE)
BERSK2_GO_Genes <- BERSK2_GO_Genes %>% rename("GO_ID.x"="BER_GO_ID", "Description.x" = "BER_Description", "GO_ID.y" = "SK2_GO_ID", "Description.y" = "SK2_Description")
BERSK2_GO_genes_surv <- merge(BERSK2_GO_Genes, survival, by="Gene")


BERSK2_GO_genes_surv_filtered <- BERSK2_GO_genes_surv %>%
  filter(!is.na(BER_GO_ID) & !is.na(SK2_GO_ID))
#JN_BER_GO_genes_filtered <- JNBERSK2_GO_genes_surv %>%
  filter(rowSums(!is.na(across(c(BER_GO_ID, JN_GO_ID, SK2_GO_ID)))) >= 2)

###making Venn of hypo/hyper methylated and up/down regulated
##mostly useless right now; need to filter a different way/think about it
BER_hyper_methy_df <- BER_merged_df_swap %>%
  filter(diff.Methy > 0)  # Assuming positive values are hyper-methylated

BER_hypo_methy_df <- BER_merged_df_swap %>%
  filter(diff.Methy < 0) 

BER_up_hyper_methy <- BER_hyper_methy_df$Gene[BER_hyper_methy_df$log2_fold_change > 0]  
BER_down_hyper_methy <- BER_hyper_methy_df$Gene[BER_hyper_methy_df$log2_fold_change < 0]
BER_only_hyper_methy <- BER_hyper_methy_df$Gene[abs(BER_hyper_methy_df$log2_fold_change) <= 0]

venn.diagram(
  x = list(BER_up_hyper_methy, BER_down_hyper_methy, BER_only_hyper_methy),
  category.names = c("Up-regulation", "Down-regulation", "Neither"),
  filename = "hyper_methylation_venn.png",
  output = TRUE ,
  imagetype="png",
  main = "Hyper-methylation",
  main.cex = 1.5,
  cex = 1.5,
  cat.cex = 1.5,
  fill = c("blue", "lightblue", "darkblue")
)



################################################################################
##doing differentially methylated LOCI now#############

setwd("C:/Users/Flinchum/Documents/Tulane Docs/Lee Lab/Results/Methyl Analysis/Dean")
## Differentially methylated LOCI
BER_minusDox_hyper_DML <- read.csv("BER_DML_peakAnno_Dox_hyper.csv")
BER_plusDox_hyper_DML <- read.csv("BER_DML_peakAnno_pDox_hyper.csv")

JN_minusDox_hyper_DML <- read.csv("JN_DML_peakAnno_Dox_hyper.csv")
JN_plusDox_hyper_DML <- read.csv("JN_DML_peakAnno_pDox_hyper.csv")

SK2_minusDox_hyper_DML <- read.csv("SK2_DML_peakAnno_Dox_hyper.csv")
SK2_plusDox_hyper_DML <- read.csv("SK2_DML_peakAnno_pDox_hyper.csv")

Tumor2911_hyper_DML <- read.csv("DML_peakAnno_primary_hyper.csv")
Tumor2912_hyper_DML <- read.csv("DML_peakAnno_recurrent_hyper.csv")

##doing gene counts for the DMLs; in two groups, the ones that are becoming hyper methylated and hypo methylated
BER_hypo_gene_counts <- BER_minusDox_hyper_DML %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n())

JN_hypo_gene_counts <- JN_minusDox_hyper_DML %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n())

SK2_hypo_gene_counts <- SK2_minusDox_hyper_DML %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n())


BER_hyper_gene_counts <- BER_plusDox_hyper_DML %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n())

JN_hyper_gene_counts <- JN_plusDox_hyper_DML %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n())

SK2_hyper_gene_counts <- SK2_plusDox_hyper_DML %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n())



BER_gene_counts <- BER_hyper_merged_DML %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n())

JN_gene_counts <- JN_hyper_merged_DML %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n())

SK2_gene_counts <- SK2_hyper_merged_DML %>% 
  group_by(SYMBOL) %>% 
  summarize(count = n())

##combining into one df with all methyl values; worked well. combined all the rows, not joined
BER_hyper_merged_DML <- bind_rows(BER_minusDox_hyper_DML %>% dplyr::select(SYMBOL, start, end, diff, pval, fdr, annotation, distanceToTSS, GENENAME),
                                  BER_plusDox_hyper_DML %>% dplyr::select(SYMBOL, start, end, diff, pval, fdr, annotation, distanceToTSS, GENENAME))

JN_hyper_merged_DML <- bind_rows(JN_minusDox_hyper_DML %>% dplyr::select(SYMBOL, start, end, diff, pval, fdr, annotation, distanceToTSS, GENENAME),
                                 JN_plusDox_hyper_DML %>% dplyr::select(SYMBOL, start, end, diff, pval, fdr, annotation, distanceToTSS, GENENAME))

SK2_hyper_merged_DML <- bind_rows(SK2_minusDox_hyper_DML %>% dplyr::select(SYMBOL, start, end, diff, pval, fdr, annotation, distanceToTSS, GENENAME),
                                  SK2_plusDox_hyper_DML %>% dplyr::select(SYMBOL, start, end, diff, pval, fdr, annotation, distanceToTSS, GENENAME))

Tumor291_hyper_merged_DML <- bind_rows(Tumor2911_hyper_DML %>% dplyr::select(SYMBOL, start, end, diff, pval, fdr, annotation, distanceToTSS, GENENAME),
                                       Tumor2912_hyper_DML %>% dplyr::select(SYMBOL, start, end, diff, pval, fdr, annotation, distanceToTSS, GENENAME))


BER_hyper_merged_DML <- BER_hyper_merged_DML %>%
  mutate(CpG_ID = str_c(start, end, sep = "-")) %>%
  relocate(CpG_ID, .after = SYMBOL) %>% # Relocate 'coordinates' to the 2nd position.
  dplyr::select(-start, -end) %>% # Remove the 'start' and 'end' columns.
  dplyr::rename(Gene = SYMBOL) %>%
  dplyr::filter(grepl("Promoter", annotation))

JN_hyper_merged_DML <- JN_hyper_merged_DML %>%
  mutate(CpG_ID = str_c(start, end, sep = "-")) %>%
  relocate(CpG_ID, .after = SYMBOL) %>% # Relocate 'coordinates' to the 2nd position.
  dplyr::select(-start, -end) %>% # Remove the 'start' and 'end' columns.
  dplyr::rename(Gene = SYMBOL) %>%
  dplyr::filter(grepl("Promoter", annotation))

SK2_hyper_merged_DML <- SK2_hyper_merged_DML %>%
  mutate(CpG_ID = str_c(start, end, sep = "-")) %>%
  relocate(CpG_ID, .after = SYMBOL) %>% # Relocate 'coordinates' to the 2nd position.
  dplyr::select(-start, -end) %>% # Remove the 'start' and 'end' columns.
  dplyr::rename(Gene = SYMBOL) %>%
  dplyr::filter(grepl("Promoter", annotation))

Tumor291_hyper_merged_DML <- Tumor291_hyper_merged_DML %>%
  mutate(CpG_ID = str_c(start, end, sep = "-")) %>%
  relocate(CpG_ID, .after = SYMBOL) %>% # Relocate 'coordinates' to the 2nd position.
  dplyr::select(-start, -end) %>% # Remove the 'start' and 'end' columns.
  dplyr::rename(Gene = SYMBOL) %>%
  dplyr::filter(grepl("Promoter", annotation))


write.csv(BER_hyper_merged_DML, "BER_hyper_merged_DML.csv")
write.csv(JN_hyper_merged_DML, "JN_hyper_merged_DML.csv")
write.csv(SK2_hyper_merged_DML, "SK2_hyper_merged_DML.csv")


##for BER
# Merge directly (no aggregation)
BER_merged_df_DML <- BER_hyper_merged_DML %>% 
  inner_join(BER_fold_change, by = "Gene")


# Create scatter plot
# Define the quadrant boundaries
log2_threshold <- 1.5
methy_threshold <- 0.2

##reversing the orientation!
BER_merged_df_DML_swap <- BER_merged_df_DML %>%
  mutate(diff = diff * -1)

# Create a new column for quadrant assignment
BER_merged_df_DML_swap <- BER_merged_df_DML_swap %>%
  mutate(Quadrant = case_when(
    log2_fold_change < -log2_threshold & diff < -methy_threshold ~ "Q1",  # Bottom left
    log2_fold_change > log2_threshold  & diff < -methy_threshold ~ "Q2",  # Top left
    log2_fold_change < -log2_threshold & diff > methy_threshold  ~ "Q3",  # bottom right
    log2_fold_change > log2_threshold  & diff > methy_threshold  ~ "Q4",  # Top right
    TRUE                                                             ~ "NS"   # Not significant
  ))


# Find top 5 points within each quadrant
BER_top_points <- BER_merged_df_DML_swap %>%
  group_by(Quadrant) %>%
  arrange(desc(abs(diff))) %>%  # Sort by absolute diff, descending
  slice_head(n = 10)  # Take top 5 points in each group

# Sort the dataframe by absolute log2 fold change and methylation difference
BER_merged_df_DML_swap <- BER_merged_df_DML_swap[order(-abs(BER_merged_df_DML_swap$log2_fold_change), -abs(BER_merged_df_DML_swap$diff)), ]

# Select the top 200 rows
BER_top_points <- head(BER_merged_df_DML_swap, 200)


ggplot(BER_merged_df_DML_swap, aes(x = diff, y = log2_fold_change, color = Quadrant)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Q1" = "blue", "Q2" = "red", "Q3" = "#006400", "Q4" = "purple", "NS" = "gray"),
                     labels = c("Q1" = "Downregulated & Hypomethylated", 
                                "Q2" = "Upregulated & Hypomethylated", 
                                "Q3" = "Downregulated & Hypermethylated",
                                "Q4" = "Upregulated & Hypermethylated",
                                "NS" = "Not Significant")) +  # Custom colors
  labs(title = "BER: Methylation Difference vs. Log2 Expression Change after EWSR1-WT1 KD",
       x = "Methylation Difference",
       y = "Log2 Expression Chang") +
  theme_minimal() + 
  geom_hline(yintercept = c(-log2_threshold, log2_threshold), linetype = "dashed", color = "gray") +  # Add horizontal lines
  geom_vline(xintercept = c(-methy_threshold, methy_threshold), linetype = "dashed", color = "gray") +   # Add vertical lines
  geom_text_repel(data = BER_top_points, aes(label = Gene), size = 3)




# Sort the dataframe by absolute log2 fold change and methylation difference
BER_merged_df_DML_swap <- BER_merged_df_DML_swap[order(-abs(BER_merged_df_DML_swap$log2_fold_change), -abs(BER_merged_df_DML_swap$diff)), ]

# Select the top 200 rows
BER_top_points <- head(BER_merged_df_DML_swap, 200)


ggplot(BER_top_points, aes(x = diff, y = log2_fold_change, color = Quadrant)) +  # Use BER_top_points here
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Q1" = "blue", "Q2" = "red", "Q3" = "#006400", "Q4" = "purple", "NS" = "gray"),
                     labels = c("Q1" = "Downregulated & Hypomethylated", 
                                "Q2" = "Upregulated & Hypomethylated", 
                                "Q3" = "Downregulated & Hypermethylated",
                                "Q4" = "Upregulated & Hypermethylated",
                                "NS" = "Not Significant")) +  # Custom colors
  labs(title = "BER: Methylation Difference vs. Log2 Expression Change after EWSR1-WT1 KD",
       x = "Methylation Difference",
       y = "Log2 Expression Chang") +
  theme_minimal() +
  geom_hline(yintercept = c(-log2_threshold, log2_threshold), linetype = "dashed", color = "gray") +  # Add horizontal lines
  geom_vline(xintercept = c(-methy_threshold, methy_threshold), linetype = "dashed", color = "gray") +  # Add vertical lines
  geom_text_repel(aes(label = Gene), size = 3)  # No need to specify data here



####trying for p-values



###trying to make venn of Gene Up/Down with Methy Up/Down; start with DMRs for now
##also not filtering for any ranges yet, will cut it down or pick thresholds after and try again
BER_gene_up <- BER_fold_change %>% filter(log2_fold_change > 1)
BER_gene_down <- BER_fold_change %>% filter(log2_fold_change < -1)

BER_gene_test <- merge(BER_gene_up, BER_gene_down, by="Gene")

BER_gene_down <- BER_gene_down %>%
  filter(Gene!= "SPATA13" & Gene!= "Metazoa_SRP" & Gene!= "Y_RNA" & Gene!= "uc_338" & Gene!= "U3" & Gene!= "SNORA81")

##reversing the orientation!
BER_merged_swap <- BER_hyper_merged %>%
  mutate(diff.Methy = diff.Methy * -1)

BER_methy_up <- BER_merged_swap %>% filter(diff.Methy > 0.05)
BER_methy_down <- BER_merged_swap %>% filter(diff.Methy < -0.05)

BER_methy_down <- na.omit(BER_methy_down)


my_colors <- brewer.pal(4, "Set1")


venn.diagram(
  x = list(BER_gene_up$Gene, BER_gene_down$Gene, BER_methy_up$Gene, BER_methy_down$Gene),
  category.names = c("Gene-Up", "Gene-Down", "Methylation-Up", "Methylation-Down"),
  filename = "BER_DEGs_DMRs.png",
  output = FALSE ,
  imagetype="png",
  main = "BER",
  main.cex = 1.5,
  cex = 1.5,
  cat.cex = 1.5,
  fill = my_colors,
  cat.pos = c(-7, 7, 0, 0),  # Adjust these angles
  cat.dist = c(0.22, 0.22, 0.11, 0.11)
)


BER_gene_down_methy_up <- merge(BER_gene_down, BER_methy_up, by="Gene")

BER_gene_down_methy_down <- merge(BER_gene_down, BER_methy_down, by="Gene")

BER_unique_to_methyUP_gene_down <- setdiff(BER_gene_down_methy_up$Gene, BER_gene_down_methy_down$Gene)

# Create a new dataframe with only the unique genes
BER_methyUP_gene_down_only <- BER_gene_down_methy_up[BER_gene_down_methy_up$Gene %in% BER_unique_to_methyUP_gene_down, ]


##JN
JN_gene_up <- JN_fold_change %>% filter(log2_fold_change > 1)
JN_gene_down <- JN_fold_change %>% filter(log2_fold_change < -1)

JN_gene_test <- merge(JN_gene_up, JN_gene_down, by="Gene")

#JN_gene_down <- JN_gene_down %>%
 # filter(Gene!= "SPATA13" & Gene!= "Metazoa_SRP" & Gene!= "Y_RNA" & Gene!= "uc_338" & Gene!= "U3" & Gene!= #"SNORA81")

##reversing the orientation!
JN_merged_swap <- JN_hyper_merged %>%
  mutate(diff.Methy = diff.Methy * -1)

JN_methy_up <- JN_merged_swap %>% filter(diff.Methy > 0.05)
JN_methy_down <- JN_merged_swap %>% filter(diff.Methy < -0.05)

JN_methy_down <- na.omit(JN_methy_down)

rows_with_na <- which(rowSums(is.na(JN_methy_up)) > 0)
print(rows_with_na)

venn.diagram(
  x = list(JN_gene_up$Gene, JN_gene_down$Gene, JN_methy_up$Gene, JN_methy_down$Gene),
  category.names = c("Gene-Up", "Gene-Down", "Methylation-Up", "Methylation-Down"),
  filename = "JN_DEGs_DMRs.png",
  output = FALSE ,
  imagetype="png",
  main = "JN",
  main.cex = 1.5,
  cex = 1.5,
  cat.cex = 1.5,
  fill = my_colors,
  cat.pos = c(-7, 7, 0, 0),  # Adjust these angles
  cat.dist = c(0.22, 0.22, 0.11, 0.11)
)

JN_gene_down_methy_up <- merge(JN_gene_down, JN_methy_up, by="Gene")

JN_gene_down_methy_down <- merge(JN_gene_down, JN_methy_down, by="Gene")

JN_unique_to_methyUP_gene_down <- setdiff(JN_gene_down_methy_up$Gene, JN_gene_down_methy_down$Gene)

# Create a new dataframe with only the unique genes
JN_methyUP_gene_down_only <- JN_gene_down_methy_up[JN_gene_down_methy_up$Gene %in% JN_unique_to_methyUP_gene_down, ]


#SK2
SK2_gene_up <- SK2_fold_change %>% filter(log2_fold_change > 1)
SK2_gene_down <- SK2_fold_change %>% filter(log2_fold_change < -1)

SK2_gene_test <- merge(SK2_gene_up, SK2_gene_down, by="Gene")

#SK2_gene_down <- SK2_gene_down %>%
# filter(Gene!= "SPATA13" & Gene!= "Metazoa_SRP" & Gene!= "Y_RNA" & Gene!= "uc_338" & Gene!= "U3" & Gene!= #"SNORA81")

##reversing the orientation!
SK2_merged_swap <- SK2_hyper_merged %>%
  mutate(diff.Methy = diff.Methy * -1)

SK2_methy_up <- SK2_merged_swap %>% filter(diff.Methy > 0.05)
SK2_methy_down <- SK2_merged_swap %>% filter(diff.Methy < -0.05)

SK2_methy_down <- na.omit(SK2_methy_down)

rows_with_na <- which(rowSums(is.na(SK2_methy_down)) > 0)
print(rows_with_na)

venn.diagram(
  x = list(SK2_gene_up$Gene, SK2_gene_down$Gene, SK2_methy_up$Gene, SK2_methy_down$Gene),
  category.names = c("Gene-Up", "Gene-Down", "Methylation-Up", "Methylation-Down"),
  filename = "SK2_DEGs_DMRs.png",
  output = FALSE ,
  imagetype="png",
  main = "SK2",
  main.cex = 1.5,
  cex = 1.5,
  cat.cex = 1.5,
  fill = my_colors,
  cat.pos = c(-7, 7, 0, 0),  # Adjust these angles
  cat.dist = c(0.22, 0.22, 0.11, 0.11)
)

SK2_gene_down_methy_up <- merge(SK2_gene_down, SK2_methy_up, by="Gene")

SK2_gene_down_methy_down <- merge(SK2_gene_down, SK2_methy_down, by="Gene")

SK2_unique_to_methyUP_gene_down <- setdiff(SK2_gene_down_methy_up$Gene, SK2_gene_down_methy_down$Gene)

# Create a new dataframe with only the unique genes
SK2_methyUP_gene_down_only <- SK2_gene_down_methy_up[SK2_gene_down_methy_up$Gene %in% SK2_unique_to_methyUP_gene_down, ]

SK2_num_unique_genes <- SK2_gene_down_methy_up %>%
  distinct(Gene) %>%
  nrow()

# Print the number of unique genes
print(SK2_num_unique_genes)



#BER DMR counts/gene
BER_DMR_counts <- BER_hyper_merged %>%
  group_by(Gene) %>%
  summarize(num_loci = n())


ggplot(BER_DMR_counts, aes(x = num_loci, fill = num_loci, group = num_loci)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(
    x = "Number of Differentially Methylated Regions",
    y = "Number of Genes",
    title = "BER: Distribution of Differentially Methylated Regions per Gene"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(BER_DMR_counts$num_loci), by = 2)) +
  scale_fill_gradient(low = "red", high = "green", trans = "sqrt")  # Add 'trans = "sqrt"'

#BER DML counts/gene
BER_DML_counts <- BER_hyper_merged_DML %>%
  group_by(Gene) %>%
  summarize(num_loci = n())

ggplot(BER_DML_counts, aes(x = num_loci, fill = num_loci, group = num_loci)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(
    x = "Number of Differentially Methylated Loci",
    y = "Number of Genes",
    title = "BER: Distribution of Differentially Methylated Loci per Gene"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(BER_DML_counts$num_loci), by = 2)) +
  scale_fill_gradient(low = "red", high = "green", trans = "log")  # Add 'trans = "sqrt"'



#JN DMR counts/gene
JN_DMR_counts <- JN_hyper_merged %>%
  group_by(Gene) %>%
  summarize(num_loci = n())


ggplot(JN_DMR_counts, aes(x = num_loci, fill = num_loci, group = num_loci)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(
    x = "Number of Differentially Methylated Regions",
    y = "Number of Genes",
    title = "JN: Distribution of Differentially Methylated Regions per Gene"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(JN_DMR_counts$num_loci), by = 2)) +
  scale_fill_gradient(low = "red", high = "green", trans = "sqrt")  # Add 'trans = "sqrt"


#SK2 DMR counts/gene
SK2_DMR_counts <- SK2_hyper_merged %>%
  group_by(Gene) %>%
  summarize(num_loci = n())


ggplot(SK2_DMR_counts, aes(x = num_loci, fill = num_loci, group = num_loci)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(
    x = "Number of Differentially Methylated Regions",
    y = "Number of Genes",
    title = "SK2: Distribution of Differentially Methylated Regions per Gene"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(SK2_DMR_counts$num_loci), by = 2)) +
  scale_fill_gradient(low = "red", high = "green", trans = "sqrt")  # Add 'trans = "sqrt"




###methy and enhancer regions
#DML

