
ft <- fread("/Users/victoriadeleray/Desktop/MPRINT_FINAL/feature_table_export/feature-table.tsv", skip = 1)

depths <- colSums(ft[, -1, with = FALSE])
depths_df <- data.frame(Sample = names(depths), SequencingDepth = depths)


median_depth <- median(depths)
depths_df <- rbind(depths_df, data.frame(Sample = "Median", SequencingDepth = median_depth))

# write to tsv
#fwrite(depths_df, "/Users/victoriadeleray/Desktop/MPRINT_FINAL/sequencing_depth_summary.tsv", sep = "\t", quote = FALSE)

FEATURE_QZA="/Users/victoriadeleray/Desktop/MPRINT_FINAL/216873_feature-table.qza"   # <-- must be .qza, not .qzv
TAXONOMY_QZA="/Users/victoriadeleray/Desktop/MPRINT_FINAL/11_3_25_taxonomy.qza"
OUTDIR="/Users/victoriadeleray/Desktop/MPRINT_FINAL/filtered_tables_11_3_25"
mkdir -p "$OUTDIR"

qiime feature-table filter-samples \
--i-table "$FEATURE_QZA" \
--p-min-frequency 10000 \
--o-filtered-table "$OUTDIR/table_min10k.qza"

qiime tools export \
--input-path "$OUTDIR/table_min10k.qza" \
--output-path "$OUTDIR/export_min10k"

NSAMP=$(biom summarize-table -i "$OUTDIR/export_min10k/feature-table.biom" | grep -m1 "Num samples:" | awk '{print $3}')
MIN_SAMPLES=$(( (NSAMP + 9) / 10 ))
echo "Samples after depth filter: $NSAMP  ->  10% prevalence threshold = $MIN_SAMPLES samples"

qiime feature-table filter-features \
--i-table "$OUTDIR/table_min10k.qza" \
--p-min-samples $MIN_SAMPLES \
--o-filtered-table "$OUTDIR/table_min10k_preval10pct.qza"

qiime tools export \
--input-path "$OUTDIR/table_min10k_preval10pct.qza" \
--output-path "$OUTDIR/export_table1"
biom convert \
-i "$OUTDIR/export_table1/feature-table.biom" \
-o "$OUTDIR/table1_min10k_preval10pct.tsv" \
--to-tsv

qiime feature-table filter-features \
--i-table "$OUTDIR/table_min10k_preval10pct.qza" \
--m-metadata-file "$TAXONOMY_QZA" \
--o-filtered-table "$OUTDIR/table_min10k_preval10pct_taxmatched.qza"

qiime taxa collapse \
--i-table "$OUTDIR/table_min10k_preval10pct_taxmatched.qza" \
--i-taxonomy "$TAXONOMY_QZA" \
--p-level 4 \
--o-collapsed-table "$OUTDIR/table_min10k_preval10pct_ORDER.qza"

qiime tools export \
--input-path "$OUTDIR/table_min10k_preval10pct_ORDER.qza" \
--output-path "$OUTDIR/export_table2"
biom convert \
-i "$OUTDIR/export_table2/feature-table.biom" \
-o "$OUTDIR/table2_min10k_preval10pct_ORDER.tsv" \
--to-tsv

qiime feature-table rarefy \
--i-table "$OUTDIR/table_min10k_preval10pct.qza" \
--p-sampling-depth 10000 \
--o-rarefied-table "$OUTDIR/table_min10k_preval10pct_rarefy10k.qza"

qiime tools export \
--input-path "$OUTDIR/table_min10k_preval10pct_rarefy10k.qza" \
--output-path "$OUTDIR/export_table3"
biom convert \
-i "$OUTDIR/export_table3/feature-table.biom" \
-o "$OUTDIR/table3_min10k_preval10pct_rarefy10k.tsv" \
--to-tsv

qiime feature-table filter-features \
--i-table "$OUTDIR/table_min10k_preval10pct_rarefy10k.qza" \
--m-metadata-file "$TAXONOMY_QZA" \
--o-filtered-table "$OUTDIR/table_min10k_preval10pct_rarefy10k_taxmatched.qza"

qiime taxa collapse \
--i-table "$OUTDIR/table_min10k_preval10pct_rarefy10k_taxmatched.qza" \
--i-taxonomy "$TAXONOMY_QZA" \
--p-level 4 \
--o-collapsed-table "$OUTDIR/table_min10k_preval10pct_rarefy10k_ORDER.qza"

qiime tools export \
--input-path "$OUTDIR/table_min10k_preval10pct_rarefy10k_ORDER.qza" \
--output-path "$OUTDIR/export_table4"
biom convert \
-i "$OUTDIR/export_table4/feature-table.biom" \
-o "$OUTDIR/table4_min10k_preval10pct_rarefy10k_ORDER.tsv" \
--to-tsv

OUTDIR="/Users/victoriadeleray/Desktop/MPRINT_FINAL/filtered_tables_11_3_25"
ALPHA_DIR="$OUTDIR/alpha_rarefy10k"
mkdir -p "$ALPHA_DIR"

qiime diversity alpha \
--i-table "$OUTDIR/table_min10k_preval10pct_rarefy10k.qza" \
--p-metric shannon \
--o-alpha-diversity "$ALPHA_DIR/shannon_rarefy10k.qza"

qiime tools export \
--input-path "$ALPHA_DIR/shannon_rarefy10k.qza" \
--output-path "$ALPHA_DIR/exported"

############################################################
alpha_tsv <- "/Users/victoriadeleray/Desktop/MPRINT_FINAL/filtered_tables_11_3_25/alpha_rarefy10k/exported/alpha-diversity.tsv"
metadata_path <- "/Users/victoriadeleray/Desktop/MPRINT_FINAL/microbiome_metadata.txt"

shannon <- fread(alpha_tsv)
metadata <- fread(metadata_path)

library(lme4)
library(lmerTest)
library(emmeans)

dat <- shannon %>%
  rename(sample_name = V1) %>%
  inner_join(metadata, by = "sample_name") %>%
  mutate(
    mouse_id = str_replace(sample_name, "D\\d+$", ""),
    across(c(time, age, ab, vaccine), \(x) factor(x))
  )

dat_mom <- dat %>% filter(age == "mom") %>% droplevels()

dat_pup <- dat %>%
  filter(age %in% c("pup", "infant")) %>% 
  droplevels()

dat_mom

model1 <- dat_mom %>%
  dplyr::mutate(ab = factor(ab, levels = c("Mock", "Amp", "Augmentin"))) %>%
  lmerTest::lmer(formula = shannon_entropy~ab+time+(1|mouse_id))
anova(model1)
summary(model1)

dat_pup

model1 <- dat_pup %>%
  dplyr::mutate(ab = factor(ab, levels = c("Mock", "Amp", "Augmentin"))) %>%
  lmerTest::lmer(formula = shannon_entropy~ab+time+(1|mouse_id))
anova(model1)
summary(model1)

dat_pup_pcv20 <- dat_pup %>%
  dplyr::filter(vaccine == "PCV 20")

model1 <- dat_pup_pcv20 %>%
  dplyr::mutate(ab = factor(ab, levels = c("Mock", "Amp", "Augmentin"))) %>%
  lmerTest::lmer(formula = shannon_entropy~ab+time+(1|mouse_id))
anova(model1)
summary(model1)



##################
feature_table <- read_tsv("/Users/victoriadeleray/Desktop/MPRINT_FINAL/filtered_tables_11_3_25/table1_min10k_preval10pct.tsv", skip = 1)

data_clean <- feature_table %>%
  column_to_rownames("#OTU ID") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("filename")

data_clean2 <- data_clean %>%
  mutate(across(where(is.numeric), as.numeric)) %>%
  column_to_rownames("filename") %>%
  dplyr::select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
  rownames_to_column("filename")
data_clean2 <- data_clean2 %>% 
  column_to_rownames("filename")

data_clr <- decostand(data_clean2, method = "rclr")
rownames(data_clr) <- rownames(data_clean2)
# PCA
data_clr_whole <- mixOmics::pca(data_clr, ncomp = 3, center = TRUE, scale = TRUE)

metadata <- read_tsv("/Users/victoriadeleray/Desktop/MPRINT_FINAL/microbiome_metadata.txt")

data_clr_whole_scores <- data.frame(data_clr_whole$variates$X) %>%
  rownames_to_column("filename") %>%
  dplyr::left_join(metadata, by = c("filename" = "sample_name"))

i <- "ab"

data_clr_whole_scores_filtered <- data_clr_whole_scores %>%
  filter(!str_detect(filename, "6mix|Blank|Pool")) %>%
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
#ggsave("//Users/victoriadeleray/Desktop/MPRINT_FINAL/filtered_tables_11_3_25/PCA_microbiome_infant_11_5_25.pdf", plot = data_clr_PCA, width = 8, height = 6, dpi = 300)


common_filenames <- intersect(rownames(data_clr), data_clr_whole_scores_filtered$filename)
data_clr2 <- data_clr[common_filenames, , drop = FALSE]
data_clr_whole_scores_filtered2 <- data_clr_whole_scores_filtered %>%
  dplyr::filter(filename %in% common_filenames)

dim(data_clr2)
dim(data_clr_whole_scores_filtered2)

dist_metabolites <- vegdist(data_clr2, method = "euclidean")
disper_donor <- betadisper(dist_metabolites, data_clr_whole_scores_filtered2$ab)
anova(disper_donor)
permutest(disper_donor, permutations = 999)

set.seed(1)
permanova <- adonis2(
  data_clr2 ~ ab * time,
  data = data_clr_whole_scores_filtered2,
  method = "euclidean",
  by = "terms"
)
permanova

feature_table <- read_tsv("/Users/victoriadeleray/Desktop/MPRINT_FINAL/filtered_tables_11_3_25/table1_min10k_preval10pct.tsv", skip = 1)

data_clean <- feature_table %>%
  column_to_rownames("#OTU ID") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("filename")

data_clean2 <- data_clean %>%
  mutate(across(where(is.numeric), as.numeric)) %>%
  column_to_rownames("filename") %>%
  dplyr::select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
  rownames_to_column("filename")
data_clean2 <- data_clean2 %>% 
  column_to_rownames("filename")

data_clr <- decostand(data_clean2, method = "rclr")
rownames(data_clr) <- rownames(data_clean2)
# PCA
data_clr_whole <- mixOmics::pca(data_clr, ncomp = 3, center = TRUE, scale = TRUE)

metadata <- read_tsv("/Users/victoriadeleray/Desktop/MPRINT_FINAL/microbiome_metadata.txt")

data_clr_whole_scores <- data.frame(data_clr_whole$variates$X) %>%
  rownames_to_column("filename") %>%
  dplyr::left_join(metadata, by = c("filename" = "sample_name"))

i <- "ab"

data_clr_whole_scores_filtered <- data_clr_whole_scores %>%
  filter(!str_detect(filename, "6mix|Blank|Pool")) %>%
  filter(age == "mom") %>%
  filter(!is.na(!!sym(i)))

data_clr_PCA <- ggscatter(data_clr_whole_scores_filtered, 
                          x = "PC1", y = "PC2", 
                          color = i, 
                          shape = "time", 
                          size = 8, alpha = 0.8,
                          title = "PCA - Rarefied rclr Microbiome (mom Only)",
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
#ggsave("//Users/victoriadeleray/Desktop/MPRINT_FINAL/filtered_tables_11_3_25/PCA_microbiome_mom_11_5_25.pdf", plot = data_clr_PCA, width = 8, height = 6, dpi = 300)

common_filenames <- intersect(rownames(data_clr), data_clr_whole_scores_filtered$filename)
data_clr2 <- data_clr[common_filenames, , drop = FALSE]
data_clr_whole_scores_filtered2 <- data_clr_whole_scores_filtered %>%
  dplyr::filter(filename %in% common_filenames)

dim(data_clr2)
dim(data_clr_whole_scores_filtered2)

dist_metabolites <- vegdist(data_clr2, method = "euclidean")
disper_donor <- betadisper(dist_metabolites, data_clr_whole_scores_filtered2$ab)
anova(disper_donor)
permutest(disper_donor, permutations = 999)

set.seed(1)
permanova <- adonis2(
  data_clr2 ~ ab * time,
  data = data_clr_whole_scores_filtered2,
  method = "euclidean",
  by = "terms"
)
permanova

#####

setwd("/Users/victoriadeleray/Desktop/MPRINT_FINAL")
feature_table <- fread("/Users/victoriadeleray/Desktop/MPRINT_FINAL/filtered_tables_11_3_25/table2_min10k_preval10pct_ORDER.tsv") %>%
  column_to_rownames("OTUID") %>%      
  as.data.frame()                      

metadata <- fread("microbiome_metadata.txt") %>%
  as.data.frame()

#metadata_filtered <- metadata %>% filter(`age-ab` %in% c("mom-Amp", "mom-Mock") & time %in% c("t=1", "t=2"))

#metadata_filtered <- metadata %>% filter(`age-ab` %in% c("mom-Amp", "mom-Mock") & time %in% c("t=3", "t=4", "t=5"))

#metadata_filtered <- metadata %>% filter(`age-ab` %in% c("mom-Augmentin", "mom-Mock") & time %in% c("t=1", "t=2"))

#metadata_filtered <- metadata %>% filter(`age-ab` %in% c("mom-Augmentin", "mom-Mock") & time %in% c("t=3", "t=4", "t=5"))

#metadata_filtered <- metadata %>% filter(`age-ab` %in% c("infant-Augmentin", "infant-Mock") & time %in% c("t=3") & vaccine == "PCV 20")

#metadata_filtered <- metadata %>% filter(`age-ab` %in% c("infant-Augmentin", "infant-Mock") & time %in% c("t=4", "t=5") & vaccine == "PCV 20")

#metadata_filtered <- metadata %>% filter(`age-ab` %in% c("infant-Amp", "infant-Mock") & time %in% c("t=3") & vaccine == "PCV 20")

metadata_filtered <- metadata %>% filter(`age-ab` %in% c("infant-Amp", "infant-Mock") & time %in% c("t=4", "t=5") & vaccine == "PCV 20")


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
ggsave("aldex_infant_amp_t45.pdf",plot = aldex, width = 10, height = 10, dpi = 300)
write.table(diverging_df, file="daldex_infant_amp_t45.tsv", sep="\t", row.names=FALSE, quote=FALSE)

############################################

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


feature_table <- fread("/Users/victoriadeleray/Desktop/MPRINT_FINAL/filtered_tables_11_3_25/table2_min10k_preval10pct_ORDER.tsv") %>%
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
diablo_perf <- perf(final_diablo, validation = "loo", dist = "centroids.dist", progressBar = TRUE, near.zero.var = TRUE, perm = 500)
diablo_perf$error.rate

plotIndiv(final_diablo, legend = TRUE, title = "DIABLO: Augmentin vs Mock")

pdf("cimDiablo_test_11_6_25.pdf", width = 12, height = 12)
cimDiablo(final_diablo, comp = 1, cutoff = 0.5)
dev.off()


circosPlot(final_diablo, cutoff = 0.7, line = FALSE, size.labels = 0.5, comp = 1,
           color.blocks = c("#4F6D7A", "#3A383F"), showIntraLinks = FALSE, legend = TRUE,
           size.variables = 0.5, size.legend = 0.5)



