##Code used to generate figures for MPRINT manuscript

setwd("/Users/victoriadeleray/Desktop/MPRINT_FINAL")

library(readxl)
library(mixOmics)
library(ggpubr)
library(vegan)
library(caret)
library(patchwork)
library(tibble)
library(tidyverse)
library(viridis)
library(pheatmap)
library(forcats)
library(ggpubr)
library(rstatix)
library(broom)
library(purrr)
library(data.table)
library(gridExtra)
library(tidyverse)
library(data.table)
library(mixOmics)
library(microbiome)
library(caret)
library(dplyr)
library(readr)
library(stringr)
library(patchwork)
library(phyloseq)
library(decontam)
library(biomformat)
library(vegan)
library(tidyverse)
library(ggplot2)
library(rlang)
library(ALDEx2)

###Clean metabolomics feature table
feature_table <- fread("mzmine_quant.csv")
feature_table <- feature_table %>%
  dplyr::rename_all(~gsub(" Peak area", "", .))
colnames(feature_table)
data_transpose <- feature_table |>
  column_to_rownames("row ID") |>    
  dplyr::select(contains(".mzML")) |> 
  t() |>  
  as.data.frame() |>  
  rownames_to_column("filename")
data_blank <- data_transpose |>  
  dplyr::filter(str_detect(pattern = "Blank", filename))
blank_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                 Mean_blank = data_blank %>% column_to_rownames("filename") %>% colMeans(),
                                 SD_blank =  data_blank %>% column_to_rownames("filename") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) %>%
  dplyr::filter(Mean_blank > 0) %>% arrange(desc(Mean_blank))
data_sample <- data_transpose |>  
  dplyr::filter(str_detect(pattern = "ME", filename)) |> 
  dplyr::filter(!str_detect(filename, "Pool"))
sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample %>% column_to_rownames("filename") %>% colMeans(),
                                  SD_sample =  data_sample %>% column_to_rownames("filename") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) %>%
  dplyr::filter(Mean_sample > 0) %>% arrange(desc(Mean_sample))
specified_features <- c("5661", "4075", "3861", "5291", "4484")
feature_to_remove <- blank_feature_info %>% 
  left_join(sample_feature_info) %>%
  dplyr::filter(Mean_blank > 0) %>%
  dplyr::mutate(Sample_Blank = Mean_sample / Mean_blank) %>%
  dplyr::filter(Sample_Blank < 5 | is.na(Sample_Blank)) %>%
  dplyr::bind_rows(blank_feature_info %>% dplyr::filter(Feature %in% specified_features)) %>%
  dplyr::distinct(Feature, .keep_all = TRUE) 
metadata <- read_csv("metabolomics_metadata.csv")
data_clean <- data_transpose |> 
  dplyr::select(-c(feature_to_remove$Feature)) |> 
  dplyr::filter(!(str_detect(filename, "6mix|Pool|Blank"))) |> 
  left_join(metadata |> dplyr::select(filename2, ATTRIBUTE_Age, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, ATTRIBUTE_time), by = c("filename" = "filename2"))
data_clean[data_clean < 0] <- 0
########################################################################################################################################################################################################################################################################################
#PCA+permanova
data_clean2 <- data_clean %>%
  select(-c("ATTRIBUTE_Age", "ATTRIBUTE_abtreatment", "ATTRIBUTE_vaccine", "ATTRIBUTE_time"))

colnames(data_clean2)

data_clr <- decostand(data_clean2 %>%  
                        dplyr::select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%  
                        column_to_rownames("filename"), method = "rclr")

data_clr_whole <- mixOmics::pca(data_clr, ncomp = 3, center = TRUE, scale = TRUE)

data_clr_whole_scores_colon <- data.frame(data_clr_whole$variates$X) %>%
  rownames_to_column("filename") %>%
  dplyr::left_join(metadata, by = c("filename" = "filename2"))
i="ATTRIBUTE_abtreatment"
data_clr_whole_scores_ab <- data_clr_whole_scores_colon %>%
  dplyr::filter(!(str_detect(filename, "6mix|Blank|Pool"))) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant") %>%  
  dplyr::filter(!is.na(!!sym(i))) 

data_clr_PCA <- data_clr_whole_scores_ab %>%
  ggscatter(x = "PC1", y = "PC2", color = i, shape = "ATTRIBUTE_time", size = 8, alpha = 0.8,
            title = "PCA - Colon",
            xlab = paste("PC1 (", round(data_clr_whole$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
            ylab = paste("PC2 (", round(data_clr_whole$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = data_clr_whole_scores_ab %>%
               group_by((!!sym(i))) %>%
               summarise_at(vars(matches("PC")), mean), 
             aes(PC1, PC2, color = (!!sym(i))), size = 7, stroke = 1, shape = 16, alpha = 1) + 
  theme(plot.title = element_text(size = 18), axis.title = element_text(size = 18),
        axis.text = element_text(size = 18)) + coord_fixed()

data_clr_PCA
#ggsave("PCA_metabolomics_infant.pdf", plot = data_clr_PCA, width = 8, height = 6, dpi = 300)
#####
common_filenames <- intersect(rownames(data_clr), data_clr_whole_scores_ab$filename)
data_clr <- data_clr[common_filenames, , drop = FALSE]
data_clr_whole_scores_ab <- data_clr_whole_scores_ab %>%
  dplyr::filter(filename %in% common_filenames)

dim(data_clr)
dim(data_clr_whole_scores_ab)

dist_metabolites <- vegdist(data_clr, method = "euclidean")
disper_donor <- betadisper(dist_metabolites, data_clr_whole_scores_ab$ATTRIBUTE_abtreatment)
anova(disper_donor)
permanova <- adonis2(dist_metabolites ~ ATTRIBUTE_vaccine*ATTRIBUTE_abtreatment*ATTRIBUTE_time, data_clr_whole_scores_ab, na.action = na.omit, by = "terms")
permanova
#####
data_clean2 <- data_clean %>%
  select(-c("ATTRIBUTE_Age", "ATTRIBUTE_abtreatment", "ATTRIBUTE_vaccine", "ATTRIBUTE_time"))

colnames(data_clean2)

data_clr <- decostand(data_clean2 %>%  
                        dplyr::select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%  
                        column_to_rownames("filename"), method = "rclr")

data_clr_whole <- mixOmics::pca(data_clr, ncomp = 3, center = TRUE, scale = TRUE)

data_clr_whole_scores_colon <- data.frame(data_clr_whole$variates$X) %>%
  rownames_to_column("filename") %>%
  dplyr::left_join(metadata, by = c("filename" = "filename2"))
i="ATTRIBUTE_abtreatment"
data_clr_whole_scores_ab <- data_clr_whole_scores_colon %>%
  dplyr::filter(!(str_detect(filename, "6mix|Blank|Pool"))) %>%
  dplyr::filter(ATTRIBUTE_Age == "mom") %>%  
  dplyr::filter(!is.na(!!sym(i))) 

data_clr_PCA <- data_clr_whole_scores_ab %>%
  ggscatter(x = "PC1", y = "PC2", color = i, shape = "ATTRIBUTE_time", size = 8, alpha = 0.8,
            title = "PCA - Colon",
            xlab = paste("PC1 (", round(data_clr_whole$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
            ylab = paste("PC2 (", round(data_clr_whole$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = data_clr_whole_scores_ab %>%
               group_by((!!sym(i))) %>%
               summarise_at(vars(matches("PC")), mean), 
             aes(PC1, PC2, color = (!!sym(i))), size = 7, stroke = 1, shape = 16, alpha = 1) + 
  theme(plot.title = element_text(size = 18), axis.title = element_text(size = 18),
        axis.text = element_text(size = 18)) + coord_fixed()

data_clr_PCA
#ggsave("PCA_metabolomics_mom.pdf", plot = data_clr_PCA, width = 8, height = 6, dpi = 300)

common_filenames <- intersect(rownames(data_clr), data_clr_whole_scores_ab$filename)
data_clr <- data_clr[common_filenames, , drop = FALSE]
data_clr_whole_scores_ab <- data_clr_whole_scores_ab %>%
  dplyr::filter(filename %in% common_filenames)

dim(data_clr)
dim(data_clr_whole_scores_ab)

dist_metabolites <- vegdist(data_clr, method = "euclidean")
disper_donor <- betadisper(dist_metabolites, data_clr_whole_scores_ab$ATTRIBUTE_abtreatment)
anova(disper_donor)
permanova <- adonis2(dist_metabolites ~ ATTRIBUTE_abtreatment*ATTRIBUTE_time, data_clr_whole_scores_ab, na.action = na.omit, by = "terms")
permanova
########################################################################################################################################################################################################################################################################################
##box plot individual feature
data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `2566`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `2566`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_vaccine == "PCV 20" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "2566"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("Tryptophanol Feature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
#ggsave("/Users/victoriadeleray/Desktop/mprint/box_plots_diablo/feature2566-tryptophanol.pdf", plot = feature_plot, width = 8, height = 6)

data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `481`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `481`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_vaccine == "PCV 20" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "481"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("Propionyl Histamine Feature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/mprint/box_plots_diablo/feature481-propionyl-histamine.pdf", plot = feature_plot, width = 8, height = 6)

data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `808`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `808`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_vaccine == "PCV 20" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "808"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("5-hydroxy tryptopholFeature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/mprint/box_plots_diablo/feature808-5-hydroxy-tryptophol", plot = feature_plot, width = 8, height = 6)

data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `805`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `805`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_vaccine == "PCV 20" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "805"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("Serotonin Feature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/mprint/box_plots_diablo/feature805-serotonin", plot = feature_plot, width = 8, height = 6)

########################################################################################################################################################################################################################################################################################
##line carnitine
values_to_keep <- c("235", "294", "658", "1466", "1754", "3129", "5916", "5451", "4937", "5913", 
                    "6508", "6255", "5305", "6467", "7117", "6809", "6539", "6258", "6528", 
                    "7601", "7328", "7009", "6804", "6808", "7705")
colnames(data_clean) <- as.character(colnames(data_clean))
data_clean_acyl <- data_clean %>% select(filename, any_of(values_to_keep))

merged_df <- data_clean_acyl %>%
  left_join(metadata, by = c("filename" = "filename2"))
filtered_df <- merged_df %>%
  filter(ATTRIBUTE_Age %in% c("mom", "infant")) %>%
  mutate(
    ATTRIBUTE_time = as.factor(ATTRIBUTE_time),
    ATTRIBUTE_abtreatment = as.factor(ATTRIBUTE_abtreatment)
  )
df_long <- filtered_df %>%
  pivot_longer(cols = all_of(values_to_keep), 
               names_to = "ScanNumber", 
               values_to = "Intensity") %>%
  mutate(Log_Intensity = log(Intensity + 1))
df_long <- df_long %>%
  mutate(ScanNumber = factor(ScanNumber, levels = values_to_keep))
df_mean <- df_long %>%
  group_by(ATTRIBUTE_time, ATTRIBUTE_abtreatment, ScanNumber, ATTRIBUTE_Age) %>%
  summarise(Mean_Log_Intensity = mean(Log_Intensity, na.rm = TRUE), .groups = "drop")
p <- ggplot(df_mean, aes(x = ATTRIBUTE_time, y = Mean_Log_Intensity, group = ScanNumber, color = ScanNumber)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  facet_grid(ATTRIBUTE_Age ~ ATTRIBUTE_abtreatment, scales = "fixed") +
  labs(title = "Mean Log-Transformed Intensity Across Timepoints",
       x = "Timepoint (ATTRIBUTE_time)",
       y = "Mean Log Intensity",
       color = "Scan Number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

print(p)
ggsave("line_graph_acyl_carnitine.pdf", plot = p, width = 20, height = 8, units = "in", dpi = 300)

##line acyl lipids
values_to_keep <- c("156", "304", "320", "427", "435", "481", "575", "769", "869", "989", "1416", "1581", "1584", "1721", "2303", "2816", "2979")

colnames(data_clean) <- as.character(colnames(data_clean))
data_clean_acyl <- data_clean %>% select(filename, any_of(values_to_keep))

merged_df <- data_clean_acyl %>%
  left_join(metadata, by = c("filename" = "filename2"))
filtered_df <- merged_df %>%
  filter(ATTRIBUTE_Age %in% c("mom", "infant")) %>%
  mutate(
    ATTRIBUTE_time = as.factor(ATTRIBUTE_time),
    ATTRIBUTE_abtreatment = as.factor(ATTRIBUTE_abtreatment)
  )
df_long <- filtered_df %>%
  pivot_longer(cols = all_of(values_to_keep), 
               names_to = "ScanNumber", 
               values_to = "Intensity") %>%
  mutate(Log_Intensity = log(Intensity + 1))
df_long <- df_long %>%
  mutate(ScanNumber = factor(ScanNumber, levels = values_to_keep))
df_mean <- df_long %>%
  group_by(ATTRIBUTE_time, ATTRIBUTE_abtreatment, ScanNumber, ATTRIBUTE_Age) %>%
  summarise(Mean_Log_Intensity = mean(Log_Intensity, na.rm = TRUE), .groups = "drop")
p <- ggplot(df_mean, aes(x = ATTRIBUTE_time, y = Mean_Log_Intensity, group = ScanNumber, color = ScanNumber)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  facet_grid(ATTRIBUTE_Age ~ ATTRIBUTE_abtreatment, scales = "fixed") +
  labs(title = "Mean Log-Transformed Intensity Across Timepoints",
       x = "Timepoint (ATTRIBUTE_time)",
       y = "Mean Log Intensity",
       color = "Scan Number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

print(p)
ggsave("line_graph_acyl_lipids.pdf", plot = p, width = 20, height = 8, units = "in", dpi = 300)

########################################################################################################################################################################################################################################################################################
values_to_keep <- c("7547", "6178", "6062", "5380", "6057", "5484", "7581", "6217", "6604", 
                    "6704", "5081", "5477", "5430", "5127", "5456", "4715", "4966", "5017", 
                    "5491", "4946", "5483", "5470", "5051", "5141", "5454", "5632", 
                    "4724", "6158", "4988", "4848", "5385", "5092", "5202", 
                    "5519", "5637", "6813", "6424", "6423", "6419", "7427")

bile_list <- read_tsv("bile_list_output.tsv")

scan_numbers <- bile_list$X.Scan.

scan_numbers <- as.character(scan_numbers)

data_clean_bile <- data_clean %>% select(filename, any_of(scan_numbers))

merged_df <- data_clean_bile %>%
  left_join(metadata, by = c("filename" = "filename2"))

filtered_df <- merged_df %>%
  filter(ATTRIBUTE_Age == "infant")

filtered_df <- filtered_df %>%
  mutate(
    ATTRIBUTE_time = as.factor(ATTRIBUTE_time),
    ATTRIBUTE_abtreatment = as.factor(ATTRIBUTE_abtreatment)
  )

df_long <- filtered_df %>%
  pivot_longer(cols = all_of(scan_numbers), 
               names_to = "ScanNumber", 
               values_to = "Intensity")

df_long <- df_long %>%
  mutate(Log_Intensity = log(Intensity + 1))

df_long <- df_long %>%
  mutate(ScanNumber = as.character(ScanNumber)) %>%
  left_join(
    bile_list %>% mutate(X.Scan. = as.character(X.Scan.)), 
    by = c("ScanNumber" = "X.Scan.")
  )

ba_order <- bile_list %>%
  distinct(X.Scan., BA) %>%
  arrange(BA) %>% 
  pull(X.Scan.)

df_long <- df_long %>%
  mutate(ScanNumber = factor(ScanNumber, levels = ba_order))

mean_intensity <- df_long %>%
  group_by(ScanNumber) %>%
  summarise(mean_intensity = mean(Log_Intensity, na.rm = TRUE))

df_long <- df_long %>%
  left_join(mean_intensity, by = "ScanNumber") %>%
  mutate(Intensity_Deviation = Log_Intensity - mean_intensity)

df_long <- df_long %>%
  complete(ATTRIBUTE_time, ATTRIBUTE_abtreatment, ScanNumber, fill = list(Intensity_Deviation = NA))

df_wide <- df_long %>%
  pivot_wider(names_from = c(ATTRIBUTE_time, ATTRIBUTE_abtreatment), 
              values_from = Intensity_Deviation) %>%
  replace_na(list(Intensity_Deviation = 0)) 

hc_scans <- hclust(dist(df_wide[, -1])) 
scan_order <- df_wide$ScanNumber[hc_scans$order]

p <- ggplot(df_long, aes(x = interaction(ATTRIBUTE_time, ATTRIBUTE_abtreatment, sep = "_"), 
                         y = ScanNumber, fill = Intensity_Deviation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Log-Transformed Scan Intensities (Clustered by BA Values)",
       x = "ATTRIBUTE_time and ATTRIBUTE_abtreatment",
       y = "Scan Number",
       fill = "Intensity\nDeviation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

print(p)
ggsave("heatmap_BA_infant.pdf", plot = p, width = 10, height = 8, units = "in", dpi = 300)

values_to_keep <- c("7547", "6178", "6062", "5380", "6057", "5484", "7581", "6217", "6604", 
                    "6704", "5081", "5477", "5430", "5127", "5456", "4715", "4966", "5017", 
                    "5491", "4946", "5483", "5470", "5051", "5141", "5454", "5632", 
                    "4724", "6158", "4988", "4848", "5385", "5092", "5202", 
                    "5519", "5637", "6813", "6424", "6423", "6419", "7427")

bile_list <- read_tsv("bile_list_output.tsv")

scan_numbers <- bile_list$X.Scan.

scan_numbers <- as.character(scan_numbers)

data_clean_bile <- data_clean %>% select(filename, any_of(scan_numbers))

merged_df <- data_clean_bile %>%
  left_join(metadata, by = c("filename" = "filename2"))

filtered_df <- merged_df %>%
  filter(ATTRIBUTE_Age == "mom")

filtered_df <- filtered_df %>%
  mutate(
    ATTRIBUTE_time = as.factor(ATTRIBUTE_time),
    ATTRIBUTE_abtreatment = as.factor(ATTRIBUTE_abtreatment)
  )

df_long <- filtered_df %>%
  pivot_longer(cols = all_of(scan_numbers), 
               names_to = "ScanNumber", 
               values_to = "Intensity")

df_long <- df_long %>%
  mutate(Log_Intensity = log(Intensity + 1))

df_long <- df_long %>%
  mutate(ScanNumber = as.character(ScanNumber)) %>%
  left_join(
    bile_list %>% mutate(X.Scan. = as.character(X.Scan.)), 
    by = c("ScanNumber" = "X.Scan.")
  )

ba_order <- bile_list %>%
  distinct(X.Scan., BA) %>%
  arrange(BA) %>% 
  pull(X.Scan.)

df_long <- df_long %>%
  mutate(ScanNumber = factor(ScanNumber, levels = ba_order))

mean_intensity <- df_long %>%
  group_by(ScanNumber) %>%
  summarise(mean_intensity = mean(Log_Intensity, na.rm = TRUE))

df_long <- df_long %>%
  left_join(mean_intensity, by = "ScanNumber") %>%
  mutate(Intensity_Deviation = Log_Intensity - mean_intensity)

df_long <- df_long %>%
  complete(ATTRIBUTE_time, ATTRIBUTE_abtreatment, ScanNumber, fill = list(Intensity_Deviation = NA))

df_wide <- df_long %>%
  pivot_wider(names_from = c(ATTRIBUTE_time, ATTRIBUTE_abtreatment), 
              values_from = Intensity_Deviation) %>%
  replace_na(list(Intensity_Deviation = 0)) 

hc_scans <- hclust(dist(df_wide[, -1])) 
scan_order <- df_wide$ScanNumber[hc_scans$order]

p <- ggplot(df_long, aes(x = interaction(ATTRIBUTE_time, ATTRIBUTE_abtreatment, sep = "_"), 
                         y = ScanNumber, fill = Intensity_Deviation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Log-Transformed Scan Intensities (Clustered by BA Values)",
       x = "ATTRIBUTE_time and ATTRIBUTE_abtreatment",
       y = "Scan Number",
       fill = "Intensity\nDeviation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

print(p)
ggsave("heatmap_BA_mom.pdf", plot = p, width = 10, height = 8, units = "in", dpi = 300)

########################################################################################################################################################################################################################################################################################

##start of microbiome analysis
########################################################################################################################################################################################################################################################################################
#shannon diversity box plots
#download tsv from alpha diversity visualization page
df <- read_tsv("infant-shannon.tsv")
df_pcv <- subset(df, vaccine == "PCV")
shannon <- ggplot(df_pcv, aes(x = ab, y = shannon_entropy, fill = ab)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # Transparent box
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +  # Add jittered points
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     comparisons = list(c("Amp", "Augmentin"),
                                        c("Amp", "Mock"),
                                        c("Augmentin", "Mock"))) +
  theme_minimal() +
  labs(title = "Shannon Diversity by AB Group (PCV only)",
       x = "Antibiotic Group", y = "Shannon Entropy") +
  theme(legend.position = "none")
shannon
df_pcv$ab <- as.factor(df_pcv$ab)
df_all <- read_tsv("/Users/victoriadeleray/Desktop/mprint/microbiome/qiita/shannon_all.tsv")
df_all$`age-ab` <- as.factor(df_all$`age-ab`)
comparisons <- combn(levels(df_all$`age-ab`), 2, simplify = FALSE)
y_limits <- range(c(df_pcv$shannon_entropy, df_all$shannon_entropy), na.rm = TRUE)
shannon <- shannon + coord_cartesian(ylim = y_limits)
shannon_all_plot <- shannon_all_plot + coord_cartesian(ylim = y_limits)
combined_plot <- shannon + shannon_all_plot +
  plot_layout(ncol = 2, guides = "collect") &
  theme(axis.title.y = element_text(angle = 90))
combined_plot
ggsave("/Users/victoriadeleray/Desktop/mprint/microbiome/qiita/plots/shannon-combined.pdf",
       plot = combined_plot, width = 10, height = 10, dpi = 300)

############################################################################################################################################################################################################################################################
qiime feature-classifier classify-sklearn  --i-classifier /Users/victoriadeleray/Desktop/mprint/microbiome/gg-13-8-99-nb-classifier.qza  --i-reads rep-seqs.qza  --o-classification taxonomy.qza

qiime feature-table rarefy \
--i-table /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/216873_feature-table.qza \
--p-sampling-depth 1100 \
--o-rarefied-table /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/216873_feature-table-rarefied-1100.qza

qiime diversity alpha \
>   --i-table /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/216873_feature-table-rarefied-1100.qza \
>   --p-metric shannon \
>   --o-alpha-diversity /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/shannon_1100.qza

qiime tools export \
>   --input-path /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/shannon_1100.qza \
>   --output-path /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/shannon_1100

shannon <- read_tsv("/Users/victoriadeleray/Desktop/MPRINT_FINAL/alpha-diversity.tsv")
metadata <- read_tsv("/Users/victoriadeleray/Desktop/MPRINT_FINAL/microbiome_metadata.txt")
metadata_filtered <- metadata %>%
  filter((grepl("^infant", `age-ab`) & vaccine == "PCV 20") | grepl("^mom", `age-ab`))
merged <- shannon %>%
  rename(sample_name = filename) %>%
  inner_join(metadata_filtered, by = "sample_name")
custom_colors <- c(
  "infant-Amp" = "#f3766e",
  "infant-Augmentin" = "#39b54a",
  "infant-Mock" = "#6f94cd",
  "mom-Amp" = "#f3766e",
  "mom-Augmentin" = "#39b54a",
  "mom-Mock" = "#6f94cd",
  "NA" = "gray"
)
comparisons <- list(
  c("infant-Amp", "infant-Augmentin"),
  c("infant-Amp", "infant-Mock"),
  c("infant-Augmentin", "infant-Mock"),
  c("mom-Amp", "mom-Augmentin"),
  c("mom-Amp", "mom-Mock"),
  c("mom-Augmentin", "mom-Mock"),
  c("infant-Amp", "mom-Amp"),
  c("infant-Augmentin", "mom-Augmentin"),
  c("infant-Mock", "mom-Mock")
)
box_shannon <- ggplot(merged, aes(x = `age-ab`, y = shannon_entropy, fill = `age-ab`)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    comparisons = comparisons,
    p.adjust.method = "fdr",
    size = 8,
    bracket.size = 1.5
  ) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("Shannon Entropy"),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
box_shannon
ggsave("/Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/shannon-entropy-rarefy-1100.pdf", plot = box_shannon, width = 10, height = 10)
box_shannon <- ggplot(merged, aes(x = `age-ab`, y = shannon_entropy, fill = `age-ab`)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    comparisons = comparisons,
    p.adjust.method = "fdr",
    size = 8,
    bracket.size = 1.5,
    y.position = c(7.6, 8.0, 8.4) # <-- adjust based on number of comparisons
  ) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Shannon Entropy",
    x = "",
    y = ""
  ) +
  coord_cartesian(ylim = c(min(merged$shannon_entropy), 8)) + # <-- tight y-axis control
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
box_shannon
ggsave("/Users/victoriadeleray/Desktop/MPRINT_FINAL/microbiome/shannon-entropy-rarefy-1100-6-26-25.pdf", plot = box_shannon, width = 10, height = 10)
############################################################################################################################################################################################################################################################

feature_table <- read_tsv("/Users/victoriadeleray/Desktop/MPRINT_FINAL/exported_rarefied/feature-table.tsv", skip = 1)

# Convert to proper format: rows = samples, columns = features
data_clean <- feature_table %>%
  column_to_rownames("#OTU ID") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("filename")

# OPTIONAL: remove low-variance features before rclr
data_clean2 <- data_clean %>%
  mutate(across(where(is.numeric), as.numeric)) %>%
  column_to_rownames("filename") %>%
  dplyr::select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
  rownames_to_column("filename")

# Apply rclr transformation
data_clr <- decostand(data_clean2 %>% column_to_rownames("filename"), method = "rclr")

# PCA
data_clr_whole <- mixOmics::pca(data_clr, ncomp = 3, center = TRUE, scale = TRUE)

metadata <- read_tsv("/Users/victoriadeleray/Desktop/MPRINT_FINAL/microbiome_metadata.txt")

data_clr_whole_scores <- data.frame(data_clr_whole$variates$X) %>%
  rownames_to_column("sample_name") %>%
  left_join(metadata, by = "sample_name")

i <- "ab"

data_clr_whole_scores_filtered <- data_clr_whole_scores %>%
  filter(!str_detect(sample_name, "6mix|Blank|Pool")) %>%
  filter(age == "infant") %>%
  filter(!is.na(!!sym(i)))

data_clr_PCA <- ggscatter(data_clr_whole_scores_filtered, 
                          x = "PC1", y = "PC2", 
                          color = i, 
                          shape = "time", 
                          size = 8, alpha = 0.8,
                          title = "PCA - Rarefied rclr Microbiome (Infant Only)",
                          xlab = paste("PC1 (", round(data_clr_whole$prop_expl_var$X[1] * 100, 1), "%)", sep = ""),
                          ylab = paste("PC2 (", round(data_clr_whole$prop_expl_var$X[2] * 100, 1), "%)", sep = ""),
                          ggtheme = theme_classic()) +
  geom_point(data = data_clr_whole_scores_filtered %>%
               group_by((!!sym(i))) %>%
               summarise(across(starts_with("PC"), mean)), 
             aes(x = PC1, y = PC2, color = (!!sym(i))), 
             size = 7, stroke = 1, shape = 16, alpha = 1) + 
  theme(plot.title = element_text(size = 18), 
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18)) +
  coord_fixed()

data_clr_PCA

# Save to file
ggsave("/Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/PCA_microbiome_infant_rarefy_1100.pdf", plot = data_clr_PCA, width = 8, height = 6, dpi = 300)


data_clr_whole_scores_filtered <- data_clr_whole_scores %>%
  filter(!str_detect(sample_name, "6mix|Blank|Pool")) %>%
  filter(age == "mom") %>%
  filter(!is.na(!!sym(i)))

data_clr_PCA <- ggscatter(data_clr_whole_scores_filtered, 
                          x = "PC1", y = "PC2", 
                          color = i, 
                          shape = "time", 
                          size = 8, alpha = 0.8,
                          title = "PCA - Rarefied rclr Microbiome (Infant Only)",
                          xlab = paste("PC1 (", round(data_clr_whole$prop_expl_var$X[1] * 100, 1), "%)", sep = ""),
                          ylab = paste("PC2 (", round(data_clr_whole$prop_expl_var$X[2] * 100, 1), "%)", sep = ""),
                          ggtheme = theme_classic()) +
  geom_point(data = data_clr_whole_scores_filtered %>%
               group_by((!!sym(i))) %>%
               summarise(across(starts_with("PC"), mean)), 
             aes(x = PC1, y = PC2, color = (!!sym(i))), 
             size = 7, stroke = 1, shape = 16, alpha = 1) + 
  theme(plot.title = element_text(size = 18), 
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18)) +
  coord_fixed()

data_clr_PCA

# Save to file
ggsave("/Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/PCA_microbiome_mom_rarefy_1100.pdf", plot = data_clr_PCA, width = 8, height = 6, dpi = 300)

############################################################################################################################################################################################################################################################
setwd("/Users/victoriadeleray/Desktop/MPRINT_FINAL")
feature_table <- fread("collapsed-feature-table.tsv") %>%
  column_to_rownames("OTUID") %>%      
  as.data.frame()                      

metadata <- fread("microbiome_metadata.txt") %>%
  as.data.frame()

metadata_filtered <- metadata %>%
  filter(`age-ab` %in% c("mom-Amp", "mom-Mock") & time %in% c("t=3", "t=4", "t=5"))
metadata_filtered <- metadata %>%
  filter(`age-ab` %in% c("mom-Augmentin", "mom-Mock") & time %in% c("t=3", "t=4", "t=5"))
metadata_filtered <- metadata %>%
  filter(`age-ab` %in% c("mom-Amp", "mom-Mock") & time %in% c("t=1", "t=2"))
metadata_filtered <- metadata %>%
  filter(`age-ab` %in% c("mom-Augmentin", "mom-Mock") & time %in% c("t=1", "t=2"))
metadata_filtered <- metadata %>%
  filter(`age-ab` %in% c("infant-Augmentin", "infant-Mock") & time %in% c("t=4", "t=5") & vaccine == "PCV 20")
metadata_filtered <- metadata %>%
  filter(`age-ab` %in% c("infant-Augmentin", "infant-Mock") & time %in% c("t=3") & vaccine == "PCV 20")
metadata_filtered <- metadata %>%
  filter(`age-ab` %in% c("infant-Amp", "infant-Mock") & time %in% c("t=4", "t=5") & vaccine == "PCV 20")
metadata_filtered <- metadata %>%
  filter(`age-ab` %in% c("infant-Amp", "infant-Mock") & time %in% c("t=3") & vaccine == "PCV 20")

intersect_names <- intersect(metadata_filtered$sample_name, colnames(feature_table))

cat("Matched samples:", length(intersect_names), "out of", nrow(metadata_filtered), "\n")

feature_table_filtered <- feature_table[, intersect_names]
metadata_filtered <- metadata_filtered %>% filter(sample_name %in% intersect_names)

conds <- metadata_filtered$`age-ab`
ncol(feature_table_filtered) == length(conds) 

aldex_res <- aldex.clr(feature_table_filtered, conds, mc.samples = 128, denom = "all", verbose = TRUE)

aldex_diff <- aldex.ttest(aldex_res, paired.test = FALSE)

aldex_summary <- aldex.effect(aldex_res)

aldex_output <- data.frame(aldex_diff, aldex_summary)

sig_features <- aldex_output %>%
  rownames_to_column("taxa") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(desc(effect))

top_taxa <- sig_features %>%
  slice_max(order_by = abs(effect), n = 100) %>%
  mutate(taxa = fct_reorder(taxa, effect))

diverging_df <- top_taxa %>%
  mutate(direction = ifelse(effect > 0, "mom-Amp", "mom-Mock"),
         effect_signed = effect,
         taxa = fct_reorder(taxa, effect_signed))

aldex <- ggplot(diverging_df, aes(x = taxa, y = effect_signed, fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("mom-Amp" = "green", "mom-Mock" = "blue")) +
  labs(title = "Differential Abundance: mom-Amp vs mom-Mock",
       x = "Taxa",
       y = "Effect Size (ALDEx2)",
       fill = "Upregulated In") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6)) 

aldex
#ggsave("aldex.pdf",plot = aldex, width = 10, height = 10, dpi = 300)
#write.table(diverging_df, file="diverging_df_infant_amp_mock_t3.tsv", sep="\t", row.names=FALSE, quote=FALSE)

############################################################################################################################################################################################################################################################
#Diablo table
library(mixOmics)
library(data.table)
library(dplyr)
library(caret)


gnps_data <- fread("/Users/victoriadeleray/Desktop/MPRINT_FINAL/merged_results_with_gnps.tsv")
scans_to_keep <- as.character(gnps_data$Scan)
colnames_data_clean <- colnames(data_clean)

columns_to_keep <- c("filename", "ATTRIBUTE_Age", "ATTRIBUTE_vaccine", intersect(colnames_data_clean, scans_to_keep))
data_clean_annotations <- data_clean[, columns_to_keep]

filtered_data_clean <- data_clean_annotations %>%
  filter(ATTRIBUTE_Age == "infant", ATTRIBUTE_vaccine == "PCV 20") %>%
  select(-starts_with("ATTRIBUTE_"))


filtered_data_clean$filename <- gsub(".mzML$", "", filtered_data_clean$filename)
filtered_data_clean$filename <- sub(".*?_", "", filtered_data_clean$filename)
filtered_data_clean$filename <- sub(".*?_", "", filtered_data_clean$filename)
filtered_data_clean <- filtered_data_clean %>% rename(sample_name = filename)


feature_table <- fread("collapsed-feature-table.tsv") %>%
  column_to_rownames("OTUID") %>%
  as.data.frame()


metadata <- fread("microbiome_metadata.txt") %>% as.data.frame()

metadata_filtered <- metadata %>%
  filter(ATTRIBUTE_all %in% c("infant-Augmentin-t=4", "infant-Augmentin-t=5", 
                              "infant-Mock-t=4", "infant-Mock-t=5"))

samples_to_keep <- metadata_filtered$sample_name

feature_table_filtered <- feature_table[, colnames(feature_table) %in% samples_to_keep, drop = FALSE]
colnames(feature_table_filtered) <- gsub("15920\\.", "", gsub("O", "E", colnames(feature_table_filtered)))


feature_table_t <- t(feature_table_filtered) %>% as.data.frame()
feature_table_t$sample_name <- rownames(feature_table_t)

common_samples <- intersect(feature_table_t$sample_name, filtered_data_clean$sample_name)

feature_table_t_sub <- feature_table_t %>%
  filter(sample_name %in% common_samples) %>%
  arrange(sample_name)

filtered_data_clean_sub <- filtered_data_clean %>%
  filter(sample_name %in% common_samples) %>%
  arrange(sample_name) %>%
  distinct(sample_name, .keep_all = TRUE)

stopifnot(identical(feature_table_t_sub$sample_name, filtered_data_clean_sub$sample_name))


X1 <- feature_table_t_sub %>%
  select(-sample_name) %>%
  mutate_all(as.numeric) %>%
  as.matrix()

X2 <- filtered_data_clean_sub %>%
  select(-sample_name) %>%
  mutate_all(as.numeric) %>%
  as.matrix()

rownames(X1) <- feature_table_t_sub$sample_name
rownames(X2) <- filtered_data_clean_sub$sample_name


grouping_df <- metadata_filtered %>%
  mutate(sample_name = gsub("15920\\.", "", gsub("O", "E", sample_name))) %>%
  filter(sample_name %in% common_samples) %>%
  arrange(sample_name) %>%
  mutate(group = case_when(
    grepl("Augmentin", ATTRIBUTE_all) ~ "Augmentin",
    grepl("Mock", ATTRIBUTE_all) ~ "Mock",
    TRUE ~ NA_character_
  ))

stopifnot(!any(is.na(grouping_df$group)))
Y <- factor(grouping_df$group)


nzv_X1 <- nearZeroVar(X1)
X1_filtered <- if (length(nzv_X1) > 0) X1[, -nzv_X1, drop = FALSE] else X1

nzv_X2 <- nearZeroVar(X2)
X2_filtered <- if (length(nzv_X2) > 0) X2[, -nzv_X2, drop = FALSE] else X2

X1_filtered <- X1_filtered[, apply(X1_filtered, 2, var) != 0, drop = FALSE]
X2_filtered <- X2_filtered[, apply(X2_filtered, 2, var) != 0, drop = FALSE]

data_list <- list(feature_table = X1_filtered, filtered_data = X2_filtered)


test_keepX <- list(
  feature_table = seq(10, 20, 2),
  filtered_data = seq(10, 20, 2)
)


design <- matrix(c(0, 1,
                   1, 0), 
                 ncol = 2, 
                 dimnames = list(names(data_list), names(data_list)))


set.seed(123)
tuned_diablo <- tune.block.splsda(X = data_list, Y = Y, ncomp = 2, 
                                  test.keepX = test_keepX, design = design, 
                                  validation = "loo", dist = "centroids.dist", near.zero.var = TRUE)


final_diablo <- block.splsda(X = data_list, Y = Y, ncomp = 2,
                             keepX = tuned_diablo$choice.keepX, design = design, near.zero.var = TRUE)


set.seed(123)
diablo_perf <- perf(final_diablo, validation = "loo", dist = "centroids.dist", progressBar = TRUE, near.zero.var = TRUE)
diablo_perf$error.rate

plotIndiv(final_diablo, legend = TRUE, title = "DIABLO: Augmentin vs Mock")

pdf("cimDiablo_test.pdf", width = 12, height = 12)
cimDiablo(final_diablo, comp = 1, cutoff = 0.5)
dev.off()


circosPlot(final_diablo, cutoff = 0.7, line = FALSE, size.labels = 0.5, comp = 1,
           color.blocks = c("#4F6D7A", "#3A383F"), showIntraLinks = FALSE, legend = TRUE,
           size.variables = 0.5, size.legend = 0.5)
############################################################################################################################################################################################################################################################
#spearman
combined_df <- cbind(
  X1_final[, c(
    "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales",
    "k__Bacteria;p__Actinobacteria;c__Coriobacteriia;o__Coriobacteriales"
  )],
  X2_final
)

cor_matrix <- cor(combined_df, method = "spearman")

cor_sub <- cor_matrix[1:2, -(1:2)]

keep_features <- colnames(cor_sub)[apply(cor_sub, 2, function(x) any(abs(x) > 0.75))]

cor_filtered <- cor_sub[, keep_features]

gnps_scans <- as.character(gnps_data$Scan)

cor_filtered_subset <- cor_filtered[, colnames(cor_filtered) %in% gnps_scans, drop = FALSE]

pdf("annotated_spearman.pdf", width = 70, height = 6)  # adjust size as needed
pheatmap(cor_filtered_subset,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         main = "Spearman Correlation (GNPS-Matched Features)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101))
dev.off()

#####
#supplementary box plots which are box plots but do are not ONLY PCV20 mice

data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `481`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `481`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "481"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("Propionyl Histamine Feature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/MPRINT_FINAL/9_12_25_plots_sup/feature481-propionyl-histamine.pdf", plot = feature_plot, width = 8, height = 6)

data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `808`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `808`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "808"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("5-hydroxy tryptopholFeature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/MPRINT_FINAL/9_12_25_plots_sup/feature808-5-hydroxy-tryptophol.pdf", plot = feature_plot, width = 8, height = 6)

data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `805`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `805`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "805"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("Serotonin Feature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/MPRINT_FINAL/9_12_25_plots_sup/feature805-serotonin", plot = feature_plot, width = 8, height = 6)

data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `2566`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `2566`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "2566"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("Tryptophanol Feature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/MPRINT_FINAL/9_12_25_plots_sup/feature2566-tryptophanol.pdf", plot = feature_plot, width = 8, height = 6)

#####
#supplementary box plots which are box plots but do are only NO PCV20


data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `481`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `481`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_vaccine == "PBS" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "481"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("Propionyl Histamine Feature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/MPRINT_FINAL/9_12_25_plots_sup/feature481-propionyl-histamine-PBS.pdf", plot = feature_plot, width = 8, height = 6)

data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `808`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `808`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_vaccine == "PBS" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "808"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("5-hydroxy tryptopholFeature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/MPRINT_FINAL/9_12_25_plots_sup/feature808-5-hydroxy-tryptophol-PBS.pdf", plot = feature_plot, width = 8, height = 6)

data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `805`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `805`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "805"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("Serotonin Feature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/MPRINT_FINAL/9_12_25_plots_sup/feature805-serotonin.pdf", plot = feature_plot, width = 8, height = 6)

data_clean2 <- data_clean %>% 
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `2566`)
data_clean3 <- data_clean2 %>%
  dplyr::select(filename, ATTRIBUTE_Age, ATTRIBUTE_time, ATTRIBUTE_abtreatment, ATTRIBUTE_vaccine, `2566`) %>%
  dplyr::filter(ATTRIBUTE_Age == "infant" & ATTRIBUTE_vaccine == "PBS" & ATTRIBUTE_time %in% c("t=4", "t=5", "t=3"))
column_name <- "2566"
custom_colors <- c("Amp" = "#f3766e", "Augmentin" = "#39b54a", "Mock" = "#6f94cd")
feature_plot <- ggplot(data_clean3, aes(x = ATTRIBUTE_abtreatment, y = .data[[column_name]], fill = ATTRIBUTE_abtreatment)) +
  geom_boxplot(width = 0.55, alpha = 0.5, lwd = 1.5, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, size = 5) +
  stat_compare_means(aes(group = ATTRIBUTE_treatment), 
                     method = "wilcox.test", 
                     label = "p.format", 
                     comparisons = list(
                       c("Amp", "Augmentin"), 
                       c("Amp", "Mock"), 
                       c("Augmentin", "Mock")
                     ),
                     p.adjust.method = "fdr",
                     size = 8,
                     bracket.size = 1.5) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = paste("Tryptophanol Feature", column_name),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 2) +
  theme(
    plot.title = element_text(size = 40, hjust = 0.5),
    axis.title = element_text(size = 24),
    axis.text.y = element_text(color = "black", size = 26),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
feature_plot
ggsave("/Users/victoriadeleray/Desktop/MPRINT_FINAL/9_12_25_plots_sup/feature2566-tryptophanol-PBS.pdf", plot = feature_plot, width = 8, height = 6)

#####
#shannon_entropy original for supplementary
df <- read_tsv("infant-shannon.tsv")
df_pcv <- df
shannon <- ggplot(df_pcv, aes(x = ab, y = shannon_entropy, fill = ab)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +  
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     comparisons = list(c("Amp", "Augmentin"),
                                        c("Amp", "Mock"),
                                        c("Augmentin", "Mock"))) +
  theme_minimal() +
  labs(title = "Shannon Diversity by AB Group (PCV only)",
       x = "Antibiotic Group", y = "Shannon Entropy") +
  theme(legend.position = "none")
shannon
ggsave("/Users/victoriadeleray/Desktop/MPRINT_FINAL/9_12_25_plots_sup/entropy-infants-all.pdf", plot = shannon, width = 8, height = 6)


################################################################################

qiime feature-classifier classify-sklearn  --i-classifier /Users/victoriadeleray/Desktop/mprint/microbiome/gg-13-8-99-nb-classifier.qza  --i-reads rep-seqs.qza  --o-classification taxonomy.qza

qiime feature-table rarefy \
--i-table /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/216873_feature-table.qza \
--p-sampling-depth 1100 \
--o-rarefied-table /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/216873_feature-table-rarefied-1100.qza

qiime diversity alpha \
>   --i-table /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/216873_feature-table-rarefied-1100.qza \
>   --p-metric shannon \
>   --o-alpha-diversity /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/shannon_1100.qza

qiime tools export \
>   --input-path /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/shannon_1100.qza \
>   --output-path /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/shannon_1100

##then make the box plot in R

qiime tools export \
>   --input-path /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/216873_feature-table-rarefied-1100.qza \
>   --output-path /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/exported_rarefied

biom convert \
>   -i /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/exported_rarefied/feature-table.biom \
>   -o /Users/victoriadeleray/Desktop/MPRINT_publish/microbiome/exported_rarefied/feature-table.tsv \
>   --to-tsv

###
qiime taxa collapse 
--i-table /Users/victoriadeleray/Desktop/mprint/microbiome/qiita/no-rarefy/216873_feature-table.qza 
--i-taxonomy /Users/victoriadeleray/Desktop/mprint/microbiome/qiita/greengenes/taxonomy.qza 
--p-level 4 
--o-collapsed-table /Users/victoriadeleray/Desktop/mprint/microbiome/qiita/collapsed-feature-table-no-rarefy-order.qza

qiime tools export  
--input-path /Users/victoriadeleray/Desktop/mprint/microbiome/qiita/no-rarefy/collapsed-feature-table-no-rarefy-order.qza  
--output-path /Users/victoriadeleray/Desktop/mprint/microbiome/qiita/no-rarefy/exported-collapsed-table

biom convert \
>   -i /Users/victoriadeleray/Desktop/mprint/microbiome/qiita/no-rarefy/exported-collapsed-table/feature-table.biom \
>   -o /Users/victoriadeleray/Desktop/mprint/microbiome/qiita/no-rarefy/exported-collapsed-table/collapsed-feature-table.tsv \
>   --to-tsv

################################################################################
taxonomy <- read_tsv("/Users/victoriadeleray/Desktop/MPRINT_FINAL/exported-taxonomy/taxonomy.tsv")

taxonomy_order <- taxonomy %>%
  mutate(Order = str_extract(Taxon, "o__[^;]+")) %>%
  select(`Feature ID`, Order)


collapsed <- read_tsv("/Users/victoriadeleray/Desktop/MPRINT_FINAL/collapsed-feature-table.tsv")

linked <- taxonomy_order %>%
  inner_join(collapsed, by = c("Order" = "OTUID"))

