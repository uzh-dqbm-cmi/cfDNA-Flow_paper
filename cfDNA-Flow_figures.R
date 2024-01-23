library(ggplot2)
library(ggfortify)
library(ggbiplot)
library(readr)
library(reshape2)
library(corrplot)
library(RColorBrewer)
library(tidyverse)
library(ggh4x)
library(pROC)
library(rstatix)

###### Figure 2A ######
output_fig2a <- "../results/Fig2A.pdf"
fig2a_df <- read_csv("../data/fig2a_df.csv")
pdf(file =  output_fig2a,
    width = 7,
    height = 7)
ggplot(fig2a_df, 
       aes(x=feature, y=variance/abs(avg), color = feature)) +
  geom_boxplot(size=1.4) +
  theme_classic() +
  scale_color_manual(values = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E7298A", "#E6AB02")) +
  labs(color = "Feature") +
  theme(text=element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x = element_blank()) +
  labs(x = "Feature", y = "Dispersion")
dev.off()

###### Figure 2B ######
output_fig2b <- "../results/Fig2B.pdf"

file_paths <- c(
  "../data/fragment_lengths.csv",
  "../data/ichorCNA.csv",
  "../data/tMAD.csv",
  "../data/normalized_coverage.csv",
  "../data/liquorice_hemato.csv",
  "../data/liquorice_lung.csv")

process_csv <- function(file_path, type) {
  data <- read_csv(file_path)
  data$settings_all <- paste(data$ref_genome, data$settings, sep="_")
  data <- data[, c("sampleID", "settings_all", "value")]
  data_wide <- spread(data, settings_all, value)
  data_wide$sampleID <- paste(data_wide$sampleID, type, sep="_")
  rownames(data_wide) <- data_wide$sampleID 
  return(data_wide)
}

types <- c("lengths", "ichorCNA", "tMAD", "normalized_coverage", "liquorice_hemato", "liquorice_lung")
processed_data_list <- list()
for (i in seq_along(file_paths)) {
  processed_data_list[[i]] <- process_csv(file_paths[i], types[i])
}
all <- do.call(rbind, processed_data_list)
all <- t(all[,-1])

# Rownames edit
current_rownames <- rownames(all)
new_rownames <- gsub("mem", "", current_rownames)
new_rownames <- gsub("True", "T", new_rownames)
new_rownames <- gsub("False", "F", new_rownames)
new_rownames <- gsub("DI26830XA", "_Strict", new_rownames)
new_rownames <- gsub("40", "_Lenient", new_rownames)
rownames(all) <- new_rownames

pca_res <- prcomp(all, scale. = TRUE)

trimming <- sub("^([03]{1,2}).*", "\\1", rownames(all))
trimming <- c(rep("Untrimmed", 16), rep("Trimmed", 16))
genome <- sub("^[03]{1,2}([^_]+)_.+", "\\1", rownames(all))
gc_bias <- sub(".*_([TF])_[^_]+$", "\\1",  rownames(all))
filtering <- sub(".*_([A-Za-z]+)$", "\\1", rownames(all))

both <- paste(filtering, gc_bias, sep="_")
pca_labels <- ifelse(grepl("T", rownames(all)), "", " ")
gc_bias <- ifelse(grepl("T", gc_bias), "GC bias corr", "GC bias uncorr")

pdf(file = output_fig2b, width = 5, height = 5)
ggbiplot(pca_res, var.axes = FALSE, labels = "") +
  geom_point(aes(color = genome, shape = both, size = trimming),
             na.rm = TRUE) + 
  scale_shape_manual(values = c(1, 0, 16, 15)) +
  scale_color_manual(values = c("#D53E4F", "#FBB4AE", "#01665E", "#66C2A5")) +
  theme_classic() +
  theme(text=element_text(size=12))
dev.off()


##### Figure 2C ######
output_directory <- "../results/"

file_paths <- c(
  "../data/fragment_lengths.csv",
  "../data/ichorCNA.csv",
  "../data/tMAD.csv",
  "../data/normalized_coverage.csv",
  "../data/liquorice_hemato.csv",
  "../data/liquorice_lung.csv"
)

generate_correlation_plots <- function(file_path, output_directory) {
  data <- read_csv(file_path)
  data$preprocessing <- paste(data$ref_genome, data$settings, sep="_")
  
  data_healthy <- data %>% filter(phenotype == "healthy_control") %>% select(c('sampleID','value','preprocessing'))
  matrix_healthy <- dcast(data_healthy, sampleID ~ preprocessing, value.var="value")
  cor_matrix_healthy <- cor(matrix_healthy[,2:length(matrix_healthy)], method = "pearson")
  
  data_cancer <- data %>% filter(phenotype == "lung_cancer") %>% select(c('sampleID','value','preprocessing'))
  matrix_cancer <- dcast(data_cancer, sampleID ~ preprocessing, value.var="value")
  cor_matrix_cancer <- cor(matrix_cancer[,2:length(matrix_cancer)], method = "pearson")
  
  base_name <- tools::file_path_sans_ext(basename(file_path))
  output_file <- paste0(output_directory, "Fig2C_", base_name, ".pdf")
  pdf(file = output_file, width = 10, height = 6) 
  
  diag(cor_matrix_healthy) <- NA
  corrplot(
    cor_matrix_healthy,
    #title = paste("\n", base_name),
    is.corr = FALSE,
    method = "circle",
    tl.cex = 0.7,
    tl.col = "black",
    type = 'lower',
    col = brewer.pal(n = 9, name = "Greens"),
    tl.srt = 45,
    #tl.pos = 'lt',
    na.label = " ",
    #cl.pos = 'n',
    tl.pos = 'n'
  )
  
  diag(cor_matrix_cancer) <- NA
  corrplot(
    cor_matrix_cancer,
    is.corr = FALSE,
    method = "circle",
    tl.cex = 0.7,
    tl.col = "black",
    type = 'upper',
    col = brewer.pal(n = 9, name = "Purples"),
    add = TRUE,
    tl.srt = 0,
    tl.pos = 'l',
    na.label = " ",
    #cl.pos = 'n',
    #tl.pos = 'n'
  )
  
  dev.off() 
}

for (input_file in file_paths) {
  generate_correlation_plots(input_file, output_directory)
}

##### Figure 3A ######
fig3a_df <- read_csv("../data/fig3a_df.csv")
pdf(file =  "../results/Fig3A.pdf",
    width = 12.5,
    height = 5)
ggplot(fig3a_df, aes(x = row_number, y = -log10(value), color = feature, group = feature, shape = Filtering)) +
  geom_point(size = 5) +
  scale_shape_manual(values = c(1, 16)) +
  geom_line(fig3a_df[fig3a_df$Filtering == 'Lenient',], mapping = aes(x = row_number, y = -log10(value),  group = feature, color = feature), alpha=0.5, size = 0.1) +
  geom_line(fig3a_df[fig3a_df$Filtering == 'Strict',], mapping = aes(x = row_number, y = -log10(value),  group = feature, color = feature), alpha=1, size = 0.1) +
  labs(x = "Preprocessing setting", y = "-log10(p)") +
  theme_minimal() +
  scale_color_manual(values = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E7298A",  "#E6AB02")) +
  scale_fill_manual(values = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E7298A", "#E6AB02")) +
  labs(fill = "Feature", col = "Feature") +
  theme(text=element_text(size=20), axis.text.x=element_blank()) 
dev.off()

##### Figure 3B #####
auc_df_plot <- read_csv("../data/auc_plot.csv")
pdf(file =  "../results/Fig3B.pdf",
    width = 5,
    height = 5)
ggplot(auc_df_plot, aes(x = feature, y = AUC, fill = feature, group = feature)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E7298A",  "#E6AB02")) +
  theme_classic() +
  theme(axis.text.x=element_blank())
dev.off()


##### Supplementary Figure 1 ######

counts_summary <- read_csv("../data/counts_summary.csv")
counts_summary$step_verbose <- factor(counts_summary$step_verbose, levels = c("Fastq", 
                                                                              "Trimming", "Alignment",
                                                                              "Filtering",
                                                                              "Final"))

pdf(file = "../results/SFig1A.pdf",
    width = 10,
    height = 5)
ggplot(counts_summary, aes(step_verbose, percent_remaining, group=ref_genome)) +   
  geom_line(aes(color=trimming), position = position_dodge(width = 0.15), 
            alpha=0.8, 
            size=0.5) +
  geom_point(aes(color=trimming), size = 2,  alpha = 0.8) + 
  scale_color_brewer(palette="Dark2") +
  theme_bw() +
  ylab("Read percentage remaining (%)") + 
  xlab("Bioinformatic preprocessing steps") +
  ggtitle("") +
  facet_nested(gc_bias ~ genome + settings_verbose) + 
  theme(panel.spacing=unit(0.1,"lines"), strip.background=element_rect(color="grey30", fill="grey90"),
        panel.border=element_rect(color="grey90"),
        axis.ticks.x=element_blank()) +
  guides(color=guide_legend(title=NULL)) +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = 45))
dev.off()

average_fragment_length <- read_csv("../data/average_fragment_length.csv")
pdf(file = "../results/SFig1B.pdf",
    width = 6,
    height = 5)
ggplot(average_fragment_length, aes(x = reads_verbose, y = avg, fill = phenotype)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.2, position = position_dodge(0.9)) +
  labs(x = "", y = "Average fragment length") +
  scale_fill_manual(values = c("healthy_control" = "#FFFF99", "lung_cancer" = "#666666")) +
  theme_classic() +
  theme(legend.position = "top") +
  guides(fill=guide_legend(title=NULL)) +
  scale_fill_brewer(palette = "Accent", labels = c('Healthy', 'Cancer')) 
dev.off()

df_summary <- read_csv("../data/discarded_reads.csv")
pdf(file = "../results/SFig1C.pdf",
    width = 9,
    height = 5)
ggplot(df_summary, aes(x = reads_verbose, y = avg, fill = phenotype)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.2, position = position_dodge(0.9)) +
  labs(x = "", y = "Average tumor fraction") +
  #scale_fill_manual(values = c("healthy_control" = "#FFFF99", "lung_cancer" = "#666666")) +
  theme_classic() +
  theme(legend.position = "top") +
  facet_nested(~size_verbose + settings_verbose) +
  guides(fill=guide_legend(title=NULL)) +
  theme(text=element_text(size=18)) +
  scale_fill_brewer(palette = "Accent", labels = c('Healthy', 'Cancer')) 
dev.off()

runtimes_cumsum <- read_csv("../data/runtimes.csv")
runtimes_cumsum$rule <- factor(runtimes_cumsum$rule, levels = c("trim", 
                                                                "bwaMem", "filterBams",
                                                                "CorrectGCBias",
                                                                "sortByCoord"))
###### Prepare Supplementary Figures ########
file_names <- c("../data/normalized_coverage.csv", 
                "../data/liquorice_hemato.csv",
                "../data/liquorice_lung.csv",
                "../data/fragment_lengths.csv", 
                "../data/ichorCNA.csv",
                "../data/tMAD.csv")
dataset_list <- list()
for (file in file_names) {
  dataset <- read.csv(file, stringsAsFactors = FALSE)
  val.col <- sub(".csv", "", sub("../data/", "", file))
  names(dataset)[names(dataset) == "value"] <- val.col
  dataset_list[[val.col]] <- subset(dataset, select=c("sampleID", "ref_genome",
                                                      "settings", "phenotype",
                                                      eval(val.col)))
}

common_column <- c("sampleID", "ref_genome", "settings", "phenotype")
merged_data <- Reduce(function(x, y) merge(x, y, by = common_column), dataset_list)

merged_data <- merged_data[!grepl("DI1230XA", merged_data$setting),]
merged_data <- separate(merged_data, col = ref_genome, into = c("trim", "build"), sep = "mem")
merged_data$gc[grepl("False", merged_data$setting)] <- FALSE
merged_data$gc[grepl("True", merged_data$setting)] <- TRUE
merged_data$filter[grepl("40", merged_data$setting)] <- "lenient"
merged_data$filter[grepl("DI26830XA", merged_data$setting)] <- "strict"
merged_data$preprocessing <- paste(merged_data$trim, merged_data$build, merged_data$settings,  sep="_")
ichor_stat.test <- merged_data %>%
  group_by(trim, build, gc, filter) %>%
  t_test(ichorCNA ~ phenotype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
ichor_stat.test$statistic <- abs(ichor_stat.test$statistic)

ichor_ttest <- ichor_stat.test %>%
  t_test(statistic ~ trim) %>%
  add_significance()
ichor_ttest$.y. <- "ichor"

tmad_stat.test <- merged_data %>%
  group_by(trim, build, gc, filter) %>%
  t_test(tMAD ~ phenotype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
tmad_stat.test$statistic <- abs(tmad_stat.test$statistic)

tmad_ttest <- tmad_stat.test %>%
  t_test(statistic ~ trim) %>%
  add_significance()
tmad_ttest$.y. <- "tmad"

length_stat.test <- merged_data %>%
  group_by(trim, build, gc, filter) %>%
  t_test(fragment_lengths ~ phenotype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
length_stat.test$statistic <- abs(length_stat.test$statistic)

length_ttest <- length_stat.test %>%
  t_test(statistic ~ trim) %>%
  add_significance()
length_ttest$.y. <- "length"

cov_stat.test <- merged_data %>%
  group_by(trim, build, gc, filter) %>%
  t_test(normalized_coverage ~ phenotype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
cov_stat.test$statistic <- abs(cov_stat.test$statistic)

cov_ttest <- cov_stat.test %>%
  t_test(statistic ~ trim) %>%
  add_significance()
cov_ttest$.y. <- "coverage"

lung_stat.test <- merged_data %>%
  group_by(trim, build, gc, filter) %>%
  t_test(liquorice_lung ~ phenotype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
lung_stat.test$statistic <- abs(lung_stat.test$statistic)

lung_ttest <- lung_stat.test %>%
  t_test(statistic ~ trim) %>%
  add_significance()
lung_ttest$.y. <- "lung"

hemato_stat.test <- merged_data %>%
  group_by(trim, build, gc, filter) %>%
  t_test(liquorice_hemato ~ phenotype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
hemato_stat.test$statistic <- abs(hemato_stat.test$statistic)

hemato_ttest <- hemato_stat.test %>%
  t_test(statistic ~ trim) %>%
  add_significance()
hemato_ttest$.y. <- "hemato"

trim_ttest <- rbind(lung_ttest, hemato_ttest, cov_ttest, length_ttest, tmad_ttest, ichor_ttest)

ichor_ttest <- ichor_stat.test %>%
  t_test(statistic ~ gc) %>%
  add_significance()
ichor_ttest$.y. <- "ichor"

tmad_ttest <- tmad_stat.test %>%
  t_test(statistic ~ gc) %>%
  add_significance()
tmad_ttest$.y. <- "tmad"

length_ttest <- length_stat.test %>%
  t_test(statistic ~ gc) %>%
  add_significance()
length_ttest$.y. <- "length"

cov_ttest <- cov_stat.test %>%
  t_test(statistic ~ gc) %>%
  add_significance()
cov_ttest$.y. <- "coverage"

lung_ttest <- lung_stat.test %>%
  t_test(statistic ~ gc) %>%
  add_significance()
lung_ttest$.y. <- "lung"

hemato_ttest <- hemato_stat.test %>%
  t_test(statistic ~ gc) %>%
  add_significance()
hemato_ttest$.y. <- "hemato"

gc_ttest <- rbind(lung_ttest, hemato_ttest, cov_ttest, length_ttest, tmad_ttest, ichor_ttest)

ichor_ttest <- ichor_stat.test %>%
  t_test(statistic ~ filter) %>%
  add_significance()
ichor_ttest$.y. <- "ichor"

tmad_ttest <- tmad_stat.test %>%
  t_test(statistic ~ filter) %>%
  add_significance()
tmad_ttest$.y. <- "tmad"

length_ttest <- length_stat.test %>%
  t_test(statistic ~ filter) %>%
  add_significance()
length_ttest$.y. <- "length"

cov_ttest <- cov_stat.test %>%
  t_test(statistic ~ filter) %>%
  add_significance()
cov_ttest$.y. <- "coverage"

lung_ttest <- lung_stat.test %>%
  t_test(statistic ~ filter) %>%
  add_significance()
lung_ttest$.y. <- "lung"

hemato_ttest <- hemato_stat.test %>%
  t_test(statistic ~ filter) %>%
  add_significance()
hemato_ttest$.y. <- "hemato"

filter_ttest <- rbind(lung_ttest, hemato_ttest, cov_ttest, length_ttest, tmad_ttest, ichor_ttest)

ichor_ttest <- ichor_stat.test %>%
  t_test(statistic ~ build) %>%
  add_significance()
ichor_ttest$.y. <- "ichor"

tmad_ttest <- tmad_stat.test %>%
  t_test(statistic ~ build) %>%
  add_significance()
tmad_ttest$.y. <- "tmad"

length_ttest <- length_stat.test %>%
  t_test(statistic ~ build) %>%
  add_significance()
length_ttest$.y. <- "length"

cov_ttest <- cov_stat.test %>%
  t_test(statistic ~ build) %>%
  add_significance()
cov_ttest$.y. <- "coverage"

lung_ttest <- lung_stat.test %>%
  t_test(statistic ~ build) %>%
  add_significance()
lung_ttest$.y. <- "lung"

hemato_ttest <- hemato_stat.test %>%
  t_test(statistic ~ build) %>%
  add_significance()
hemato_ttest$.y. <- "hemato"

build_ttest <- rbind(lung_ttest, hemato_ttest, cov_ttest, length_ttest, tmad_ttest, ichor_ttest)

alltests <- rbind(trim_ttest[, 1:8], filter_ttest[, 1:8], gc_ttest[, 1:8], build_ttest[, 1:8]) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance()

all_stat_test <- rbind(ichor_stat.test, tmad_stat.test, length_stat.test, hemato_stat.test, lung_stat.test, cov_stat.test)

mapping_dict <- setNames(c("length", "ichorCNA", "tMAD", "hemato", "lung", "coverage"), c("fragment_lengths", "ichorCNA", "tMAD", "liquorice_hemato", "liquorice_lung", "normalized_coverage"))

# Use match to substitute values
all_stat_test$`.y.` <- mapping_dict[as.character(all_stat_test$`.y.`)]

all_stat_test$`.y.` <- factor(all_stat_test$`.y.`, levels=c("ichorCNA", "tMAD",
                                                            "length", "hemato",
                                                            "lung", "coverage"))

mapping_dict <- setNames(c("no", "yes"), c("0", "30"))

# Use match to substitute values
all_stat_test$trim <- mapping_dict[as.character(all_stat_test$trim)]

###### Supplementary Figure 2 #######
pdf(file =  "../results/SFig2.pdf",
    width = 12,
    height = 8)
ggplot(all_stat_test, aes(x = trim, y = statistic, fill = `.y.`, alpha = trim)) +
  geom_violin() +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  geom_point(position = position_jitter(width = 0.2), size = 1, shape = 21) +
  scale_fill_manual(values = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02","#E7298A","#E6AB02")) +
  facet_wrap(.~`.y.`, scales = "free",) +
  theme_classic()
dev.off()

##### Supplementary Figure 3 #######
pdf(file =  "../results/SFig3.pdf",
    width = 12,
    height = 8)
ggplot(all_stat_test, aes(x = build, y = statistic, fill = `.y.`, alpha = build)) +
  geom_violin() +
  scale_alpha_discrete(range = c(0.1, 0.9)) +
  geom_point(position = position_jitter(width = 0.2), size = 1, shape = 21) +
  scale_fill_manual(values = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02","#E7298A","#E6AB02")) +
  facet_wrap(.~`.y.`, scales = "free",) +
  theme_classic()
dev.off()

##### Supplementary Figure 4 ######
pdf(file =  "../results/SFig4.pdf",
    width = 12,
    height = 8)
ggplot(all_stat_test, aes(x = gc, y = statistic, fill = `.y.`, alpha = gc)) +
  geom_violin() +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  geom_point(position = position_jitter(width = 0.2), size = 1, shape = 21) +
  scale_fill_manual(values = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02","#E7298A","#E6AB02")) +
  facet_wrap(.~`.y.`, scales = "free",) +
  theme_classic()
dev.off()

##### Supplementary Figure 5 #####
pdf(file =  "../results/SFig5.pdf",
    width = 12,
    height = 8)
ggplot(all_stat_test, aes(x = filter, y = statistic, fill = `.y.`, alpha = filter)) +
  geom_violin() +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  geom_point(position = position_jitter(width = 0.2), size = 1, shape = 21) +
  scale_fill_manual(values = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02","#E7298A","#E6AB02")) +
  facet_wrap(.~`.y.`, scales = "free",) +
  theme_classic()
dev.off()

##### Supplementary Figure 6 ######

pdf(file =  "../results/SFig6.pdf",
    width = 10,
    height = 4)
ggplot(runtimes_cumsum, aes(x=rule, y=csum/3600, colour=coloring, group = groups)) + 
  geom_point(size=2) + 
  geom_line(size=1, alpha=0.6) + 
  theme_classic() + 
  facet_grid(~genome+trimming) +
  ylab("Hour") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylim(0, 15)
dev.off()

###########