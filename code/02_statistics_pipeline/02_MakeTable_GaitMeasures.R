## Make Stat Table for Kinematic Results
# Performs linear mixed-effects modeling and planned contrasts on gait kinematic measures.
# dv ~ group * cue + age + (1 | subject); Type III ANOVA with cluster-robust vcov (vcovCR, CR2).
#
# requires steps_lmer_kinematic_run1.xlsx from matlab pipeline scripts and club_contrast.r 
#
# Notes:
# - stride_mean is rescaled (cm/height) and speed_mean to cm before modeling.
#
# LM 050426
##

rm(list = ls())
library(dplyr)
library(lmerTest)
library(car)
library(tidyverse)
library(readxl)
library(writexl)
library(emmeans)
library(sjPlot)
library(ggeffects)
library(performance)
library(ggsignif)

# set up paths
currentDir <- "E:/clab/DoD-Gait/results" # <-- Change this to your local path

dataDir <- file.path(currentDir, "processed_data")
outputDir <- file.path(currentDir, "output")
funcDir <- file.path(currentDir,"a1_statistics_pipeline","func")

# load in functions
source(file.path(funcDir, "club_contrast.r"))
source(file.path(funcDir, "calculate_tstat.r"))

df <- read_excel(file.path(dataDir, "steps_lmer_kinematic_run1.xlsx"))
df$group <- factor(df$group)
df$subgroup <- factor(df$subgroup, levels = c("hc", "nofog", "fog"))
df$subgroup <- relevel(df$subgroup, "hc")
df$group <- relevel(df$group, "hc")
df$cue <- gsub("nocue", "No Cue", df$cue) # change label
df$cue <- gsub("auditory", "Auditory", df$cue)
df$cue <- gsub("visual", "Visual", df$cue)
df$cue <- factor(df$cue, levels = unique(df$cue))
df$cue <- relevel(df$cue, "No Cue")
df <- df %>% # rescale stride and speed to cm 
  mutate(
    stride_mean =  stride_mean * 100 / height,
    speed_mean =  speed_mean * 100
  )
## MUST SET CONSTRASTS for ANOVA
contrasts(df$cue) <- contr.sum
contrasts(df$group) <- contr.sum
contrasts(df$subgroup) <- contr.sum

## make models and contrasts ##################################################
#dv_list <- c("speed_mean", "speed_cv",
#  "stride_mean", "stride_cv",
#  "cadence_mean", "cadence_cv")

dv_list <- c("speed_mean", "stride_mean", "cadence_cv")

# function for extract stats
process_model <- function(dv, df, skip_random_effect) {

  # build formula
  if (skip_random_effect) {
    formula <- as.formula(paste(dv, "~ subgroup * cue "))
  } else {
    formula <- as.formula(paste(dv, "~ subgroup * cue + age + (1 | subject)"))
  }
  model <- lmer(formula, data = df)
  model_vcov_cr <- vcovCR(model, type = "CR2", cluster = model.frame(model)$subject)
  anova_result <- car::Anova(model, type = 3, vcov = model_vcov_cr)

  anova_df <- anova_result %>%
    as.data.frame() %>%
    rownames_to_column(var = "Predictors") %>%
    mutate(variable = dv) %>%
    select(variable, everything())

  # Planned contrasts
  cue_emm <- emmeans(model, ~ cue)
  cue_main <- club_contrast(model, model_vcov_cr, cue_emm)

  cue_pair_emm <- emmeans(model, ~ cue | subgroup)
  cue_pair <- club_contrast(model, model_vcov_cr, cue_pair_emm)

  group_emm <- emmeans(model, ~ subgroup)
  group_main <- club_contrast(model, model_vcov_cr, group_emm)

  group_pair_emm <- emmeans(model, ~ subgroup | cue)
  group_pair <- club_contrast(model, model_vcov_cr, group_pair_emm)
  
  contrast_types <- list(
    cue_main,
    cue_pair,
    group_main,
    group_pair
  )

  names(contrast_types) <- c("cue_main", "cue_pair", "group_main", "group_pair")
  
  
  contrast_results <- lapply(contrast_types, function(emm) {
    emm %>%
      as.data.frame() %>%
      mutate(
        p_val.fdr = p.adjust(p_val, method = "fdr"),
        signif.fdr = case_when(
          p_val.fdr < 0.05  ~ "*",
          TRUE ~ ""
        ),
        variable = dv
      ) %>%
      select(variable, everything())
  })

  list(
    model = model,
    anova = anova_df,
    contrasts = contrast_results,
    cue_emm = cue_emm,
    group_emm = group_emm
  )
}

results <- vector("list", length(dv_list))
names(results) <- dv_list
for (i in seq_along(dv_list)) {
  dv <- dv_list[i]
  results[[i]] <- process_model(dv, df, skip_random_effect = FALSE)
}
## end of models ##

## make ANOVA table ##
results_ANOVA_combined <- bind_rows(
  lapply(names(results), function(dv) {
    results[[dv]]$anova %>%
      mutate(
        variable = dv,
        signif = case_when(
          `Pr(>Chisq)` < 0.05  ~ "*",
          TRUE ~ ""
        )
      ) %>%
      select(variable, everything())
  })
)

results_ANOVA_combined <- results_ANOVA_combined %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
#file_path <- file.path(outputDir, "Table-Kinematic_ANOVA_all.csv")
#write.csv(results_ANOVA_combined, file_path, row.names = FALSE)

## combine EMM
results_cue_combined <- bind_rows(
  lapply(names(results), function(dv) {
    results[[dv]]$cue_emm %>%
      as.data.frame() %>%
      mutate(variable = dv)
  })
)

results_group_combined <- bind_rows(
  lapply(names(results), function(dv) {
    results[[dv]]$group_emm %>%
      as.data.frame() %>%
      mutate(variable = dv)
  })
)

## make Cue Main Effects table ##
results_cue_main_combined <- bind_rows(
  lapply(names(results), function(dv) {
    results[[dv]]$contrasts$cue_main %>%
      mutate(
        variable = dv,
        signif = case_when(
          p_val < 0.05  ~ "*",
          TRUE ~ ""
        )
      ) %>%
      select(variable, everything())
  })
)

results_cue_main_combined <- calculate_tstat(results_cue_main_combined, "Est", "SE") %>%
  mutate(across(where(is.character), ~ replace_na(., "")))

#file_path <- file.path(outputDir, "Table-Kinematic_cue_main_all.csv")
#write.csv(results_cue_main_combined, file_path, row.names = FALSE)

## make Cue Pair Contrast table ##
results_cue_pair_combined <- bind_rows(
  lapply(names(results), function(dv) {
    results[[dv]]$contrasts$cue_pair %>%
      mutate(variable = dv)
  })
) %>%
  extract(
    contrast,
    into = c("comp1", "comp2"),
    regex = "([a-zA-Z ]+) - (.+)",
    remove = FALSE
  )

## add t-statistic
results_cue_pair_combined <- calculate_tstat(results_cue_pair_combined, "Est", "SE") %>%
  mutate(across(where(is.character), ~ replace_na(., "")))

## clean up table for for supplementary table
results_cue_pair_combined_clean <- results_cue_pair_combined %>%
  select(variable, subgroup, contrast, Est, SE, t_stat, df, CI_L, CI_U, p_val, p_val.fdr)

## Flip contrasts where No Cue is the first comp
results_cue_pair_combined_clean <- results_cue_pair_combined_clean %>%
  separate(contrast, into = c("c1", "c2"), sep = " - ", remove = FALSE) %>%
  mutate(
    flip = c1 == "No Cue",
    orig_CI_L = CI_L,
    orig_CI_U = CI_U,
    new_c1 = ifelse(flip, c2, c1),
    new_c2 = ifelse(flip, c1, c2),
    new_contrast = paste(new_c1, new_c2, sep = " - "),
    Est  = ifelse(flip, -Est, Est),
    CI_L  = ifelse(flip, -orig_CI_U, orig_CI_L),
    CI_U  = ifelse(flip, -orig_CI_L, orig_CI_U),
    t_stat = ifelse(flip, -t_stat, t_stat)
  ) %>%
  select(-contrast, -c1, -c2, -new_c1, -new_c2, -flip, -orig_CI_L, -orig_CI_U) %>%
  rename(contrast = new_contrast) %>%
  relocate(contrast, .before = 2)

results_cue_pair_combined_clean <- rename(
  results_cue_pair_combined_clean,
  Variable = variable,
  Subgroup = subgroup,
  p = p_val,
  p.fdr = p_val.fdr,
  Contrast = contrast,
  CI95_L = CI_L,
  CI95_U = CI_U
)

results_cue_pair_combined_clean$Variable <- dplyr::recode(
  results_cue_pair_combined_clean$Variable,
  "speed_mean"  = "Gait Speed (cm/s)",
  "stride_mean" = "Stride Length (cm/height)",
  "cadence_cv"  = "Cadence Variability (%)"
)

results_cue_pair_combined_clean$Subgroup <- dplyr::recode(
  results_cue_pair_combined_clean$Subgroup,
  "hc"  = "HC",
  "nofog" = "PD-NF",
  "fog" = "PD-F"
)

results_cue_pair_combined_clean <- results_cue_pair_combined_clean %>%
  mutate( across(where(is.numeric), ~ round(.x, 3)))

file_path <- file.path(outputDir, "Table-Kinematics_Pairwise_Contrast.csv")
write.csv(results_cue_pair_combined_clean, file_path, row.names = FALSE)


## make significance table ################################################
### Coding for geom_signif
## HC = 1; PD-NF = 2; PD-F = 3
## Auditory = -.25; No Cue = 0; Visual = .25

results_p <- results_cue_pair_combined %>%
  select(variable, subgroup, comp1, comp2, Est, SE, t_stat, df, CI_L, CI_U, p_val, signif, p_val.fdr, signif.fdr) %>%
  filter(p_val < 0.05) %>%
  mutate(
    group_val = case_when(
      subgroup == "hc" ~ 1,
      subgroup == "nofog" ~ 2,
      subgroup == "fog" ~ 3,
      TRUE ~ NA_real_
    ),
    comp1_val = case_when(
      comp1 == "Auditory" ~ -0.25,
      comp1 == "No Cue"   ~ 0.0,
      comp1 == "Visual"   ~ 0.25,
      TRUE ~ NA_real_
    ),
    comp2_val = case_when(
      comp2 == "Auditory" ~ -0.25,
      comp2 == "No Cue"   ~ 0.0,
      comp2 == "Visual"   ~ 0.25,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    min = pmin(comp1_val + group_val, comp2_val + group_val, na.rm = TRUE),
    max = pmax(comp1_val + group_val, comp2_val + group_val, na.rm = TRUE)
  ) %>%
  select(variable, subgroup, comp1, comp2, Est, SE, t_stat, df, CI_L, CI_U, p_val, signif, p_val.fdr, signif.fdr, min, max)

file_path <- file.path(outputDir, "Table-Kinematics_4plot_Cue_Pairwise.csv")
write.csv(results_p, file_path, row.names = FALSE)

## create group contrast ##################################################
results_group_main_combined <- bind_rows(
  lapply(names(results), function(dv) {
    results[[dv]]$contrasts$group_main %>%
      mutate(
        variable = dv,
      ) %>%
      select(variable, everything())
  })
)

results_group_main_combined <- results_group_main_combined %>%
  extract(
    contrast,
    into = c("comp1", "comp2"),
    regex = "(.+) - (.+)"
  )

##% add t-statistic
results_group_main_combined <- calculate_tstat(results_group_main_combined, "Est", "SE") %>%
  mutate(across(where(is.character), ~ replace_na(., "")))

## save to excel
output_list <- list(
  ANOVA = results_ANOVA_combined,
  Cue_EMM = results_cue_combined,
  Group_EMM = results_group_combined,
  Cue_Main = results_cue_main_combined,
  Cue_Pair = results_cue_pair_combined,
  Group_Main = results_group_main_combined
)

file_path <- file.path(outputDir, "Table-Kinematics_Analysis_All.xlsx")
write_xlsx(output_list, path = file_path)

## INTERACTION ANALYSIS ##
## pairwise testing for speed model
m1 <- lmer(speed_mean ~ subgroup * cue + age + (1 | subject), df)
m1_vcov_cr <- vcovCR(m1, type = "CR2", cluster = model.frame(m1)$subject)
m1_result <- car::Anova(m1, type = 3, vcov = m1_vcov_cr)
emm_interaction1 <- emmeans(m1, ~ subgroup * cue)
contrast_results1 <- club_contrast(m1, m1_vcov_cr, emm_interaction1, inter = TRUE)
contrast_df1 <- as.data.frame(contrast_results1) %>%
  mutate(
    p_val.fdr = p.adjust(p_val, method = "fdr"),
    signif.fdr = case_when(
      p_val.fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## pairwise testing for stride model
m2 <- lmer(stride_mean ~ subgroup * cue + age + (1 | subject), df)
m2_vcov_cr <- vcovCR(m2, type = "CR2", cluster = model.frame(m2)$subject)
m2_result <- car::Anova(m2, type = 3, vcov = m2_vcov_cr)
emm_interaction2 <- emmeans(m2, ~ subgroup * cue)
contrast_results2 <- club_contrast(m2, m2_vcov_cr, emm_interaction2, inter = TRUE)
contrast_df2 <- as.data.frame(contrast_results2) %>%
  mutate(
    p_val.fdr = p.adjust(p_val, method = "fdr"),
    signif.fdr = case_when(
      p_val.fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

contrast_all <- rbind(
  contrast_df1 %>% mutate(model = "speed"),
  contrast_df2 %>% mutate(model = "stride")
)

  contrast_all <- calculate_tstat(contrast_all, "Est", "SE") %>%
    mutate(across(where(is.character), ~ replace_na(., "")))

file_path <- file.path(outputDir, "Table-Kinematics_Interaction_Contrasts.csv")
write.csv(contrast_all, file_path, row.names = FALSE)

## baseline comparison for cadence_mean check if there are group differences in cadence_mean in No Cue condition only
df_summary <- df %>%
  filter(cue == "No Cue") %>%
  group_by(subgroup) %>%
  summarise(
    mean_cadence = mean(cadence_mean, na.rm = TRUE),
    sd_cadence   = sd(cadence_mean, na.rm = TRUE),
    max_cadence = max(cadence_mean, na.rm = TRUE),
    n            = n()
  )

df %>% filter(cue == "No Cue") %>% group_by(subgroup) %>% summarise(n= n())

df_summary <- df_summary %>%
  mutate(
    mean_cadence = sprintf("%.2f", mean_cadence),
    sd_cadence   = sprintf("%.2f", sd_cadence),
    max_cadence   = sprintf("%.2f", max_cadence)
  )

df_noCue <- df %>%
  filter(cue == "No Cue")
  
anova_result <- aov(cadence_mean ~ subgroup, data = df_noCue)
summary(anova_result)

