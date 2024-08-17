# Set Up ------------------------------------------------------------------

# Edit required
getwd()
setwd("C:/Users/fishm/OneDrive/Documents/MSc EEC/Thesis/Results/all_analysis/MSc_EEC_Thesis_Tyler_Christian")

library(tidyverse)
library(factoextra)
library(cluster)
library(lme4)
library(patchwork)
library(MCMCglmm)
library(stats)
library(parallel)

# Data --------------------------------------------------------------------

# Note: all data files have been formatted for data tidying and for optimal use in later analysis.

# Extracted acoustic parameters using avisoft
best_calls <- read.csv("avisoft_data.csv") %>%
  separate(audio_marker,
           into = c("audio_marker", "before_after"),
           sep = "_(?=a|b)") %>%
  mutate(
    audio_marker = as.factor(audio_marker),
    before_after = as.factor(before_after),
    site = as.factor(site),
    call = as.factor(call),
    segment = as.factor(segment)
  )

# Amounts of calls per audio clip segmented using Koe
koe_data <- read.csv("koe_data.csv", na.strings = "") %>%
  mutate(call_type = sequence) %>%
  filter(quality != "INVALID") %>%
  group_by(filename) %>%
  separate(filename,
           into = c("audio_marker", "before_after"),
           sep = "_(?=a|b)") %>%
  separate_longer_delim(cols = call_type, delim = "-") %>%
  mutate(across(call_type, ~ str_remove_all(., '"')),
         across(call_type, ~ str_remove_all(., '_'))) %>%
  group_by(audio_marker, before_after, site, duration) %>%
  count(call_type) %>%
  pivot_wider(names_from = call_type, values_from = n) %>%
  replace(is.na(.), 0) %>%
  select(!c("NA")) %>%
  mutate(
    audio_marker = as.factor(audio_marker),
    before_after = as.factor(before_after),
    site = as.factor(site),
    site = gsub(x = site, pattern = "ascot", replacement = "mainland"),
    site = gsub(x = site, pattern = "essex", replacement = "mainland")
  ) %>%
  mutate( 
    CC_per_time = CC / duration, # standardizing call frequencies according to audio clip duration
    DW_per_time = DW / duration,
    HQ_per_time = HQ / duration,
    CA_per_time = CA / duration,
    NC_per_time = NC / duration,
    TQ_per_time = TQ / duration
  )

# All field data recorded on Lundy and the Mainland
behaviour_data <- read.csv("behavioural_data.csv") %>%
  replace(is.na(.), 0) %>%
  mutate(
    date = as.factor(date),
    site = as.factor(site),
    taxidermy = as.factor(taxidermy),
    audio_marker = as.factor(audio_marker),
    time = as.factor(time),
    taxidermy_present = as.factor(taxidermy_present),
    site = gsub(x = site, pattern = "ascot", replacement = "mainland"), # combining data from the 2 mainland sites to ensure a reasonable sample size
    site = gsub(x = site, pattern = "essex", replacement = "mainland")
  ) %>%
  filter(usable_stimuli == "x") %>% # filterning for observations marked as usable (i.e. meet the stimuli requirements)
  select(
    !c(
      trial,
      observation,
      number_male,
      number_female,
      audio_start:angle_video_finish,
      wild_pigeon:wild_other_passerine,
      usable_stimuli,
      notes
    )
  )

# Merging call frequencies with behavioural data
combined_data <- merge(behaviour_data, koe_data)

# used in analyzing the taxidermy data
combined_data_long <- combined_data %>% 
  group_by(audio_marker, before_after, site) %>% 
  mutate(wild_starling = ifelse(wild_starling == 0,
                                NA,
                                "starling"),
         wild_corvid = ifelse(wild_corvid == 0,
                              NA,
                              "corvid"),
         wild_gull = ifelse(wild_gull == 0,
                            NA,
                            "gull"),
         wild_raptor = ifelse(wild_raptor == 0,
                              NA,
                              "raptor"),
         human_disturbance = ifelse(human_disturbance == 0,
                                    NA,
                                    "human")) %>% 
  pivot_longer(cols = c(wild_starling:human_disturbance, taxidermy_present), 
               names_to = "stimuli",
               values_to = "stimuli_type") %>%
  filter(stimuli_type != "closed") %>% 
  na.omit() %>% 
  mutate(taxidermy_wild = ifelse(stimuli_type == "sparrow_taxidermy",
                                 "taxidermy",
                                 ifelse(stimuli_type == "woodpecker_taxidermy",
                                        "taxidermy",
                                        ifelse(stimuli_type == "starling_taxidermy",
                                               "taxidermy",
                                               ifelse(stimuli_type == "sparrowhawk_taxidermy",
                                                      "taxidermy",
                                                      ifelse(stimuli_type == "empty_taxidermy",
                                                             "taxidermy",
                                                             "wild"))))),
         control_threat = ifelse(stimuli_type == "empty_taxidermy",
                                  "control",
                                  ifelse(stimuli_type == "woodpecker_taxidermy",
                                         "control",
                                         "threatening"))) %>% 
  ungroup()

# Used in linear mixed models to make interpretation easier
combined_data_long_binary_ba <- combined_data_long %>% 
  mutate(before_after = gsub(x = before_after, pattern = "b", replacement = "0"),
         before_after = gsub(x = before_after, pattern = "a", replacement = "1"))

# used for creating Fig. 3.3.1
calls_long <- read.csv("koe_data.csv", na.strings = "") %>%
  mutate(call_type = sequence) %>%
  filter(quality != "INVALID") %>%
  group_by(filename) %>%
  separate(filename,
           into = c("audio_marker", "before_after"),
           sep = "_(?=a|b)") %>%
  separate_longer_delim(cols = call_type, delim = "-") %>%
  mutate(
    across(call_type, ~ str_remove_all(., '"')),
    across(call_type, ~ str_remove_all(., '_')),
    site = gsub(x = site, pattern = "ascot", replacement = "mainland"),
    site = gsub(x = site, pattern = "essex", replacement = "mainland")
  ) %>%
  filter(site == "lundy") %>% 
  select(audio_marker, before_after, site, call_type, duration) %>%
  group_by(audio_marker, before_after, site) %>%
  count(call_type) %>%
  ungroup()

# used for creating Fig. 3.3.1
calls_long_mean_no_stimuli <- merge(calls_long, behaviour_data) %>% 
  filter(taxidermy_present != c("empty_taxidermy", "woodpecker_taxidermy")) %>% 
  mutate(wild_starling = ifelse(wild_starling == 0,
                                NA,
                                "wild_starling"),
         wild_corvid = ifelse(wild_corvid == 0,
                              NA,
                              "wild_corvid"),
         wild_gull = ifelse(wild_gull == 0,
                            NA,
                            "wild_gull"),
         wild_raptor = ifelse(wild_raptor == 0,
                              NA,
                              "wild_raptor"),
         human_disturbance = ifelse(human_disturbance == 0,
                                    NA,
                                    "human_disturbance")) %>% 
  pivot_longer(cols = c(wild_starling:human_disturbance), names_to = "remove_me", values_to = "stimuli") %>% 
  select(c(site, before_after, call_type, stimuli, n)) %>% 
  na.omit() %>% 
  mutate(site = as.factor(site),
         before_after = as.factor(before_after),
         call_type = as.factor(call_type),
         stimuli = as.factor(stimuli)) %>% 
  filter(stimuli != "closed") %>% 
  group_by(before_after, call_type) %>% 
  summarise("mean" = mean(n),
            "length" = length(n),
            "sd" = sd(n),
            "se" = sd/sqrt(length)) %>% # calculating standard error for error bars
  mutate(call_type = replace_na(call_type, "none"),
         before_after = gsub(pattern = "a", replacement = "After", x = before_after),
         before_after = gsub(pattern = "b", replacement = "Before", x = before_after)) %>% 
  ungroup() %>% 
  add_row(before_after = "Before", call_type = "HQ", mean = 0) %>% 
  add_row(before_after = "Before", call_type = "TQ", mean = 0) %>% 
  replace(is.na(.), 0)

# Cluster Analysis All Data --------------------------------------------------------------

# constructing a hierarchical cluster analysis using all of the appropriate extracted acoustic parameters

# Aggregate the data, filtering for lundy and removing HQ (only 1 sample but I don't want to delete data just incase)
aggregated_data <- best_calls %>% 
  filter(site == "lundy")

# View the aggregated data
print(aggregated_data)

# Standardize the data
data_scaled <- scale(aggregated_data %>%
                       select(c(7,9:14,16:21,23:27)))

# Calculate the distance matrix
dist_matrix <- dist(data_scaled, 
                    method = "euclidean")

# Perform hierarchical clustering
hc <- hclust(dist_matrix, 
             method = "ward.D2")

# cophenetic correlation coefficient: 
coph_corr <- cor(cophenetic(hc), dist(aggregated_data))
coph_corr # acceptable

### The cophenetic correlation coefficient measures how faithfully a dendrogram preserves the pairwise distances between the original data points. A higher value indicates a better representation of the data.
# 0-0.5 bad
# 0.5-0.75 acceptable
# 0.75-1.0 very good

# Cluster Analysis Iterations --------------------------------------------------

### Warning: DO NOT RUN!!! ###
### May take 30 min to run!!!
### You have been warned!!!

# Uses parallel processing (much faster than loops but still slow) to identify the optimal acoustic parameters.

# Define the column groups
column_groups <- list(
  c(7, 14, 21),
  c(9, 16, 23),
  c(10, 17, 24),
  c(11, 18, 25),
  c(12, 19, 26),
  c(13, 20, 27)
)

# Optional column
optional_column <- 6

# Function to calculate the cophenetic correlation coefficient for a given set of columns
calculate_cophenetic_corr <- function(data, columns) {
  data_scaled <- scale(data %>% select(all_of(columns)))
  dist_matrix <- dist(data_scaled, method = "euclidean")
  hc <- hclust(dist_matrix, method = "ward.D2")
  cophenetic_corr <- cor(cophenetic(hc), dist(dist_matrix))
  return(cophenetic_corr)
}

# Generate all combinations of columns from the groups
generate_combinations <- function(groups, optional) {
  all_combinations <- expand.grid(lapply(groups, function(group) group))
  colnames(all_combinations) <- paste0("V", 1:ncol(all_combinations))
  combinations_with_optional <- rbind(
    cbind(Optional = optional, all_combinations),
    cbind(Optional = NA, all_combinations)
  )
  return(combinations_with_optional)
}

# Find the best combination of columns
best_combination <- NULL
best_cophenetic_corr <- -Inf

# Generate all possible combinations
column_combinations <- generate_combinations(column_groups, optional_column)
total_combinations <- nrow(column_combinations)
current_combination <- 0

# Function to process each combination
process_combination <- function(i) {
  combination <- na.omit(as.numeric(column_combinations[i, ]))
  cophenetic_corr <- calculate_cophenetic_corr(aggregated_data, combination)
  return(list(combination = combination, cophenetic_corr = cophenetic_corr))
}

# Counter function to print progress
counter_function <- function(i) {
  current_combination <<- current_combination + 1
  cat("Processing combination", current_combination, "of", total_combinations, "\n")
  return(process_combination(i))
}

# Run the combinations in parallel using parLapply
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

# Load necessary libraries and export objects to each node in the cluster
clusterEvalQ(cl, {
  library(dplyr)
  library(stats)
})
clusterExport(cl, list("column_combinations", "calculate_cophenetic_corr", "aggregated_data", "process_combination", "current_combination", "total_combinations"))

# Run the combinations in parallel
results <- parLapply(cl, 1:total_combinations, function(i) {
  current_combination <<- current_combination + 1
  cat("Processing combination", current_combination, "of", total_combinations, "\n")
  return(process_combination(i))
})

stopCluster(cl)

# Find the best result
for (result in results) {
  if (result$cophenetic_corr > best_cophenetic_corr) {
    best_cophenetic_corr <- result$cophenetic_corr
    best_combination <- result$combination
  }
}

# Print the best combination of columns and the corresponding cophenetic correlation coefficient
cat("Best combination of columns:", best_combination, "\n")
cat("Best cophenetic correlation coefficient:", best_cophenetic_corr, "\n")

### Output:
# Best combination of columns: 14 23 17 18 26 27 
# Best cophenetic correlation coefficient: 0.8334661 

## This is much better than using all acoustic parameters

# Optimized Cluster Analysis ----------------------------------------------------

# Aggregate the data, filtering for lundy
aggregated_data_optimised <- best_calls %>% 
  filter(site == "lundy")

# View the aggregated data
print(aggregated_data_optimised)

# Standardize the data
data_scaled_optimised <- scale(aggregated_data_optimised %>%
                                 select(c(14,23,17,18,26,27)))

# Calculate the distance matrix (using optimal methods suggested by the literature)
dist_matrix_optimised <- dist(data_scaled_optimised, 
                              method = "euclidean")

# Perform hierarchical clustering
hc_optimised <- hclust(dist_matrix_optimised, 
                       method = "ward.D2")

# Plotting the Deprogram Fig. 3.4.1 -------------------------------------------------

# Set colour groups per call
aggregated_data_optimised <- aggregated_data_optimised %>% 
  mutate(color_call = ifelse(call == "CA",
                             "#88ccee",
                             ifelse(call == "CC",
                                    "#661100",
                                    ifelse(call == "DW",
                                           "#117733",
                                           ifelse(call == "HQ",
                                                  "#999933",
                                                  ifelse(call == "NC",
                                                         "#e31a1c",
                                                         "#aa4499" ### "TQ"
                                                  )
                                           )
                                    )
                             )
  )
  ) %>% 
  mutate(color_call = as.factor(color_call))

# Enhanced dendrogram visualization with color-coded nodes by call (Fig. 3.2.1)
dendrogram_optimised <-  # Creating base dendrogram
  fviz_dend(
    hc_optimised,
    rect = FALSE,
    cex = 0.3,
    lwd = 0.5,
    k_colors = "black",
    labels_track_height = 0.8,
    show_labels = TRUE,
    label_cols = aggregated_data_optimised$color_call
  ) +
  theme(plot.title = element_blank()) +
  annotate( # adding manually interpreted groups as a visual aid
    "rect",
    xmin = 0,
    xmax = 44,
    ymin = -1.8,
    ymax = -0.9,
    alpha = .3,
    fill = "#88ccee",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 22,
    y = -1.35,
    size = 2.5,
    label = "CA"
  ) +
  annotate(
    "rect",
    xmin = 44,
    xmax = 95,
    ymin = -1.8,
    ymax = -0.9,
    alpha = .3,
    fill = "#661100",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 69.5,
    y = -1.35,
    size = 2.5,
    label = "CC"
  ) +
  annotate(
    "rect",
    xmin = 95,
    xmax = 103,
    ymin = -1.8,
    ymax = -0.9,
    alpha = .3,
    fill = "#e31a1c",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 99,
    y = -1.35,
    size = 2.5,
    label = "NC"
  ) +
  annotate(
    "rect",
    xmin = 103,
    xmax = 118,
    ymin = -1.8,
    ymax = -0.9,
    alpha = .3,
    fill = "#88ccee",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 110.5,
    y = -1.35,
    size = 2.5,
    label = "CA"
  ) +
  annotate(
    "rect",
    xmin = 118,
    xmax = 125,
    ymin = -1.8,
    ymax = -0.9,
    alpha = .3,
    fill = "#117733",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 121.5,
    y = -1.35,
    size = 2.5,
    label = "DW"
  ) +
  annotate(
    "rect",
    xmin = 125,
    xmax = 137,
    ymin = -1.8,
    ymax = -0.9,
    alpha = .3,
    fill = "#999933",
    color = "black",
    lwd = 0.3
  )  +
  annotate(
    "text",
    x = 131,
    y = -1.35,
    size = 2.5,
    label = "HQ"
  ) +
  annotate(
    "rect",
    xmin = 137,
    xmax = 142,
    ymin = -1.8,
    ymax = -0.9,
    alpha = .3,
    fill = "#aa4499",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 139.5,
    y = -1.35,
    size = 2.5,
    label = "TQ"
  ) +
  annotate(
    "rect",
    xmin = 142,
    xmax = 146,
    ymin = -1.8,
    ymax = -0.9,
    alpha = .3,
    fill = "#e31a1c",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 144,
    y = -1.35,
    size = 2.5,
    label = "NC"
  ) +
  annotate(
    "rect",
    xmin = 146,
    xmax = 159,
    ymin = -1.8,
    ymax = -0.9,
    alpha = .3,
    fill = "#aa4499",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 152.5,
    y = -1.35,
    size = 2.5,
    label = "TQ"
  ) +
  annotate(
    "rect",
    xmin = 159,
    xmax = 167,
    ymin = -1.8,
    ymax = -0.9,
    alpha = .3,
    fill = "#e31a1c",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 163,
    y = -1.35,
    size = 2.5,
    label = "NC"
  ) +
  annotate(
    "text",
    x = -5.1,
    y = -1.35,
    size = 2.5,
    label = "Group"
  )

dendrogram_optimised

ggsave("Fig_3.4.1_dendrogram_optimised.jpg", width = 10.4, height = 6, units = "in")

# Linear Mixed Models of Call Rates ----------------------------------------------------------------------

CC_model <- lmer(CC_per_time ~ before_after + (1|taxidermy_wild/control_threat) + (1|site) + total_sparrows, data = combined_data_long_binary_ba)
NC_model <- lmer(NC_per_time ~ before_after + (1|taxidermy_wild/control_threat) + (1|site) + total_sparrows, data = combined_data_long_binary_ba)
TQ_model <- lmer(TQ_per_time ~ before_after + (1|taxidermy_wild/control_threat) + (1|site) + total_sparrows, data = combined_data_long_binary_ba)
HQ_model <- lmer(HQ_per_time ~ before_after + (1|taxidermy_wild/control_threat) + (1|site) + total_sparrows, data = combined_data_long_binary_ba)
CA_model <- lmer(CA_per_time ~ before_after + (1|taxidermy_wild/control_threat) + (1|site) + total_sparrows, data = combined_data_long_binary_ba)
DW_model <- lmer(DW_per_time ~ before_after + (1|taxidermy_wild/control_threat) + (1|site) + total_sparrows, data = combined_data_long_binary_ba)

summary(CC_model)
summary(NC_model)
summary(TQ_model) # significant
summary(HQ_model)
summary(CA_model)
summary(DW_model)

# Before/After plot Fig. 3.3.1 -------------------------------------------------------

calls_long_mean_no_stimuli$call_type <-
  factor(calls_long_mean_no_stimuli$call_type,
         levels = c("CC", "NC", "TQ", "HQ", "CA", "DW"))

before_after_plot <- calls_long_mean_no_stimuli %>% 
  ggplot(aes(
    x = before_after,
    y = mean,
    group = call_type,
    color = call_type
  )) +
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.04,
                lwd = 0.2) +
  ylab("Mean Calls") +
  labs(color = "Call Type") +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(limits = rev) +
  scale_color_manual(values = c("#88ccee","#661100","#117733","#999933","#e31a1c","#aa4499"))

before_after_plot

ggsave("Fig_3.3.1_before_after_plot.jpg", height = 7.5, width = 6.8, units = "in")

# Linear Models for Taxidermy Effectiveness -------------------------------

# Filter for starling stimuli (wild+taxidermy)
combined_data_long_starling <- combined_data_long %>% 
  filter(stimuli_type == c("starling", "starling_taxidermy"),
         before_after == "a")

CC_model_stimuli <- lm(CC_per_time ~ taxidermy_wild, data = combined_data_long_starling)
NC_model_stimuli <- lm(NC_per_time ~ taxidermy_wild, data = combined_data_long_starling)
TQ_model_stimuli <- lm(TQ_per_time ~ taxidermy_wild, data = combined_data_long_starling)
HQ_model_stimuli <- lm(HQ_per_time ~ taxidermy_wild, data = combined_data_long_starling)
CA_model_stimuli <- lm(CA_per_time ~ taxidermy_wild, data = combined_data_long_starling)
DW_model_stimuli <- lm(DW_per_time ~ taxidermy_wild, data = combined_data_long_starling)

# all insignificant
summary(CC_model_stimuli)
summary(NC_model_stimuli)
summary(TQ_model_stimuli) # no data
summary(HQ_model_stimuli) 
summary(CA_model_stimuli)
summary(DW_model_stimuli)

# Linear Models of Call Rates per Site (Island Effect) ----------------------------------

CC_model_site <- lmer(CC_per_time ~ site + (1|before_after), data = combined_data_long_binary_ba)
NC_model_site <- lmer(NC_per_time ~ site + (1|before_after), data = combined_data_long_binary_ba)
TQ_model_site <- lmer(TQ_per_time ~ site + (1|before_after), data = combined_data_long_binary_ba)
HQ_model_site <- lmer(HQ_per_time ~ site + (1|before_after), data = combined_data_long_binary_ba)
CA_model_site <- lmer(CA_per_time ~ site + (1|before_after), data = combined_data_long_binary_ba)
DW_model_site <- lmer(DW_per_time ~ site + (1|before_after), data = combined_data_long_binary_ba)

summary(CC_model_site)
summary(NC_model_site) # significantly more on mainland
summary(TQ_model_site) 
summary(HQ_model_site) 
summary(CA_model_site)
summary(DW_model_site)

# Island Effect PCA set-up -------------------------------------------------------

call_pca_data <- best_calls %>% 
  mutate(site = gsub(x = site, pattern = "ascot", replacement = "Mainland"),
         site = gsub(x = site, pattern = "essex", replacement = "Mainland"),
         site = gsub(x = site, pattern = "lundy", replacement = "Lundy"))

call_pca <- princomp(call_pca_data[, c(7,9:14,16:21,23:27)], # remove duration and peak amplitude
                     cor = TRUE)

### Latent Root Criterion:

eigenvalues <- call_pca$sdev^2
eigenvalues
# 2 components

### The Scree Plot Criterion:

plot(call_pca, type="lines", ylim=c(0,10))
# 3 components

### The Relative Percent Variance Criterion:

summary(call_pca)
# 2 components represent 85.6% of variance

# therefore use the first 2 components

### Interpreting the component structure:
loadings(call_pca)

# Peak frequency (end) - Component 2
# Band width (mean) - Component 1
# 25% quartile (end) - Component 2
# 50% quartile (end) - Component 1
# 75% quartile (mean) - Component 1
# Entropy (mean) - Component 1

# Island Effect PCA Plot Fig. 3.1.1 -------------------------------------------------------

# Insufficient data for TQ, HQ and DW PCA

# Preparing PCA data
call_pca_data_ALL <- best_calls %>% 
  mutate(site = gsub(x = site, pattern = "ascot", replacement = "Mainland"),
         site = gsub(x = site, pattern = "essex", replacement = "Mainland"),
         site = gsub(x = site, pattern = "lundy", replacement = "Lundy"))

call_pca_data_CC <- best_calls %>% 
  mutate(site = gsub(x = site, pattern = "ascot", replacement = "Mainland"),
         site = gsub(x = site, pattern = "essex", replacement = "Mainland"),
         site = gsub(x = site, pattern = "lundy", replacement = "Lundy")) %>% 
  filter(call == "CC")

call_pca_data_NC <- best_calls %>% 
  mutate(site = gsub(x = site, pattern = "ascot", replacement = "Mainland"),
         site = gsub(x = site, pattern = "essex", replacement = "Mainland"),
         site = gsub(x = site, pattern = "lundy", replacement = "Lundy")) %>% 
  filter(call == "NC")

call_pca_data_CA <- best_calls %>% 
  mutate(site = gsub(x = site, pattern = "ascot", replacement = "Mainland"),
         site = gsub(x = site, pattern = "essex", replacement = "Mainland"),
         site = gsub(x = site, pattern = "lundy", replacement = "Lundy")) %>% 
  filter(call == "CA")

# Creating PCA
call_pca_ALL <- princomp(call_pca_data_ALL[, c(7,9:14,16:21,23:27)], # remove duration and peak amplitude
                     cor = TRUE)

call_pca_CC <- princomp(call_pca_data_CC[, c(7,9:14,16:21,23:27)], # remove duration and peak amplitude
                     cor = TRUE)

call_pca_NC <- princomp(call_pca_data_NC[, c(7,9:14,16:21,23:27)], # remove duration and peak amplitude
                     cor = TRUE)

call_pca_CA <- princomp(call_pca_data_CA[, c(7,9:14,16:21,23:27)], # remove duration and peak amplitude
                     cor = TRUE)

# Plotting PCA
call_pca_plot_CC <- fviz_pca_biplot(call_pca_CC, 
                                 repel = TRUE,
                                 geom.ind = "point",
                                 geom.var = "none",
                                 ellipse.level = 0.95,
                                 col.var = "black", 
                                 labelsize = 4,
                                 habillage = call_pca_data_CC$site,
                                 addEllipses = TRUE,
                                 alpha.ind = 0.5,
                                 title = "CC") +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("darkorchid3", "orange4")) +
  scale_color_manual(values = c("darkorchid3", "orange4"))

call_pca_plot_NC <- fviz_pca_biplot(call_pca_NC, 
                                    repel = TRUE,
                                    geom.ind = "point",
                                    geom.var = "none",
                                    ellipse.level = 0.95,
                                    col.var = "black", 
                                    labelsize = 4,
                                    habillage = call_pca_data_NC$site,
                                    addEllipses = TRUE,
                                    alpha.ind = 0.5,
                                    title = "NC") +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("darkorchid3", "orange4")) +
  scale_color_manual(values = c("darkorchid3", "orange4"))

call_pca_plot_CA <- fviz_pca_biplot(call_pca_CA, 
                                    repel = TRUE,
                                    geom.ind = "point",
                                    geom.var = "none",
                                    ellipse.level = 0.95,
                                    col.var = "black", 
                                    labelsize = 4,
                                    habillage = call_pca_data_CA$site,
                                    addEllipses = TRUE,
                                    alpha.ind = 0.5,
                                    title = "CA") +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("darkorchid3", "orange4")) +
  scale_color_manual(values = c("darkorchid3", "orange4"))

call_pca_plot_ALL_call <- fviz_pca_biplot(call_pca_ALL, 
                                          repel = TRUE,
                                          geom.ind = "point",
                                          geom.var = "none",
                                          ellipse.level = 0.95,
                                          col.var = "black", 
                                          labelsize = 4,
                                          habillage = call_pca_data$site,
                                          addEllipses = TRUE,
                                          alpha.ind = 0.5,
                                          title = "All Calls") +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("darkorchid3", "orange4")) +
  scale_color_manual(values = c("darkorchid3", "orange4"))

(call_pca_plot_CC + call_pca_plot_NC) / 
  (call_pca_plot_CA + call_pca_plot_ALL_call) +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

ggsave("Fig_3.1.1_PCA.jpg", height = 6.2, width = 6.8, units = "in")

