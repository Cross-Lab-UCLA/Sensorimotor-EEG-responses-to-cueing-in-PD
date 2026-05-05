## Make stat tables for the LSM and RSM side models
# Performs linear mixed-effects modeling and planned contrasts on gait kinematic measures.
# dv ~ subgroup * cue + age + (1 | subject); Type III ANOVA with cluster-robust vcov (vcovCR, CR2).
#
# requires erspData from matlab pipeline scripts and club_contrast.r
#
# LM 050426
##

# set up libraries
rm(list = ls())
library(lmerTest)
library(tidyverse)
library(readxl)
library(writexl)
library(emmeans)
library(clubSandwich)

# set up paths
currentDir <- "E:/clab/DoD-Gait/results" # <-- Change this to your local path

dataDir <- file.path(currentDir, "processed_data")
outputDir <- file.path(currentDir, "output")
funcDir <- file.path(currentDir, "a1_statistics_pipeline", "func")

# load in functions
source(file.path(funcDir, "club_contrast.r"))
source(file.path(funcDir, "calculate_tstat.r"))

# dependent variables
# dvs <- c("betaSupp_selfBase_pChange", "betaMod_selfBase_pChange", "muMod_selfBase_pChange")
dvs <- c("betaSupp_selfBase_pChange")

# filenames
RSM_filename <- "RSM_erspData.csv"
LSM_filename <- "LSM_erspData.csv"

drop_subjects <- c("")

clean_df <- function(filename, dataDir, drop_subjects) {
  df <- read_csv(file.path(dataDir, filename), show_col_types = FALSE)
  n_subjects <- df %>%
    summarise(n = n_distinct(subject)) %>%
    pull(n)
  message("Total subjects in ", filename, ": ", n_subjects)

  df <- df %>%
    filter(!subject %in% drop_subjects) %>%
    mutate(
      group = factor(group),
      subgroup = relevel(factor(subgroup, levels = c("hc", "nofog", "fog")), "hc"),
      cue = factor(case_when(
        cue == "nocue" ~ "No Cue",
        cue == "auditory" ~ "Auditory",
        cue == "visual" ~ "Visual",
        TRUE ~ cue
      ), levels = c("No Cue", "Auditory", "Visual"))
    )

  return(df)
}

df_RSM <- clean_df(RSM_filename, dataDir, drop_subjects)
df_LSM <- clean_df(LSM_filename, dataDir, drop_subjects)
output_list <- list()

# run loop
for (dv in dvs) {
  message("Running analysis for: ", dv)

  # set formula
  formula1 <- as.formula(paste(dv, "~ cue * subgroup + age + (1 | subject)"))

  options(contrasts = c("contr.sum", "contr.sum")) # ensure contrasts for type 3 ANOVA

  RSM_subgroup_counts <- df_RSM %>%
    distinct(subject, subgroup) %>%
    count(subgroup)
  for (i in seq_len(nrow(RSM_subgroup_counts))) {
    message(
      "RSM Subgroup '", RSM_subgroup_counts$subgroup[i],
      "', n = ", RSM_subgroup_counts$n[i], " subjects."
    )
  }

  RSM <- lmer(formula1, data = df_RSM)
  RSM_vcov_cr <- vcovCR(RSM, type = "CR2", cluster = model.frame(RSM)$subject)
  RSM_a <- car::Anova(RSM, type = 3, vcov = RSM_vcov_cr)

  RSM_anova_df <- RSM_a %>%
    as.data.frame() %>%
    mutate(model = paste0("RSM ", dv)) %>%
    tibble::rownames_to_column(var = "term")

  emm_cueCompByGroup <- emmeans(RSM, ~ cue | subgroup)

  RSM_cue_pairs_by_subgroup <- club_contrast(RSM, RSM_vcov_cr, emm_cueCompByGroup) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  RSM_emm_cue <- emmeans(RSM, ~cue)
  RSM_cue_main <- club_contrast(RSM, RSM_vcov_cr, RSM_emm_cue) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  RSM_emm_group <- emmeans(RSM, ~subgroup)
  RSM_group_main <- club_contrast(RSM, RSM_vcov_cr, RSM_emm_group) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  # LSM MODEL
  LSM_subgroup_counts <- df_LSM %>%
    distinct(subject, subgroup) %>% # Ensure one row per subject
    count(subgroup)
  for (i in seq_len(nrow(LSM_subgroup_counts))) {
    message(
      "LSM Subgroup '", LSM_subgroup_counts$subgroup[i],
      "', n = ", LSM_subgroup_counts$n[i], " subjects."
    )
  }

  LSM <- lmer(formula1, data = df_LSM)
  LSM_vcov_cr <- vcovCR(LSM, type = "CR2", cluster = model.frame(LSM)$subject)
  LSM_a <- car::Anova(LSM, type = 3, vcov = LSM_vcov_cr)

  LSM_anova_df <- LSM_a %>%
    as.data.frame() %>%
    mutate(model = paste0("LSM ", dv)) %>%
    tibble::rownames_to_column(var = "term")

  emm_cueCompByGroup <- emmeans(LSM, ~ cue | subgroup)
  LSM_cue_pairs_by_subgroup <- club_contrast(LSM, LSM_vcov_cr, emm_cueCompByGroup) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  LSM_emm_cue <- emmeans(LSM, ~cue)
  LSM_cue_main <- club_contrast(LSM, LSM_vcov_cr, LSM_emm_cue) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  LSM_emm_group <- emmeans(LSM, ~subgroup)
  LSM_group_main <- club_contrast(LSM, LSM_vcov_cr, LSM_emm_group) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  # combine anova tables from LSM and RSM
  anova_combined <- bind_rows(LSM_anova_df, RSM_anova_df)
  anova_combined <- anova_combined %>%
    select(model, everything()) %>%
    mutate(
      signif = case_when(
        `Pr(>Chisq)` < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  # combine emmean cue
  cue_combined <- bind_rows(
    as.data.frame(LSM_emm_cue) %>% mutate(variable = paste("LSM", dv)),
    as.data.frame(RSM_emm_cue) %>% mutate(variable = paste("RSM", dv))
  )

  group_combined <- bind_rows(
    as.data.frame(LSM_emm_group) %>% mutate(variable = paste("LSM", dv)),
    as.data.frame(RSM_emm_group) %>% mutate(variable = paste("RSM", dv))
  )

  # saveFile <- file.path(outputDir, paste0("Table-ANOVA_", dv, ".csv"))
  # write_csv(anova_combined, saveFile)

  ## planned contrasts ##
  cue_pairwise_combined <- bind_rows(
    LSM_cue_pairs_by_subgroup %>% mutate(variable = paste("LSM", dv)),
    RSM_cue_pairs_by_subgroup %>% mutate(variable = paste("RSM", dv))
  ) %>%
    mutate(across(where(is.character), ~ replace_na(., "")))

  cue_pairwise_combined <- cue_pairwise_combined %>%
    extract(
      contrast,
      into = c("comp1", "comp2"),
      regex = "([a-zA-Z ]+) - (.+)",
      remove = FALSE
    )

  # file_path <- file.path(outputDir, paste0("Table-Cue_Pairwise_", dv, ".csv"))
  # write.csv(cue_pairwise_combined, file_path, row.names = FALSE)

  # create table for contrasts that are signficant for plotting
  results_p <- cue_pairwise_combined %>%
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
        comp1 == "No Cue" ~ 0.0,
        comp1 == "Visual" ~ 0.25,
        TRUE ~ NA_real_
      ),
      comp2_val = case_when(
        comp2 == "Auditory" ~ -0.25,
        comp2 == "No Cue" ~ 0.0,
        comp2 == "Visual" ~ 0.25,
        TRUE ~ NA_real_
      )
    ) %>%
    mutate(
      min = pmin(comp1_val + group_val, comp2_val + group_val, na.rm = TRUE),
      max = pmax(comp1_val + group_val, comp2_val + group_val, na.rm = TRUE)
    ) %>%
    select(variable, subgroup, comp1, comp2, Est, SE, df, CI_L, CI_U, p_val, signif, p_val.fdr, signif.fdr, min, max)

  file_path <- file.path(outputDir, paste0("Table-", dv, "_4plot_Cue_Pairwise.csv"))
  write.csv(results_p, file_path, row.names = FALSE)

  ## main effects
  cue_main_combined <- bind_rows(
    LSM_cue_main %>% mutate(variable = paste("LSM", dv)),
    RSM_cue_main %>% mutate(variable = paste("RSM", dv))
  )

  group_main_combined <- bind_rows(
    LSM_group_main %>% mutate(variable = paste("LSM", dv)),
    RSM_group_main %>% mutate(variable = paste("RSM", dv))
  )

  LSM_subgroup_counts$side <- "LSM"
  RSM_subgroup_counts$side <- "RSM"
  merged_counts <- bind_rows(LSM_subgroup_counts, RSM_subgroup_counts)
  merged_counts <- merged_counts %>% select(side, everything())

  ### add t-stat
  cue_main_combined <- calculate_tstat(cue_main_combined,"Est", "SE") %>%
    mutate(across(where(is.character), ~ replace_na(., "")))

  group_main_combined <- calculate_tstat(group_main_combined, "Est", "SE") %>%
    mutate(across(where(is.character), ~ replace_na(., "")))

  cue_pairwise_combined <- calculate_tstat(cue_pairwise_combined, "Est", "SE") %>%
    mutate(across(where(is.character), ~ replace_na(., "")))
  
  ### save to excel
  output_list[[dv]] <- list(
    ANOVA = anova_combined,
    Cue_Est = cue_combined,
    Group_Est = group_combined,
    Cue_Main = cue_main_combined,
    Cue_Pair = cue_pairwise_combined,
    Group_Main = group_main_combined,
    Subject_n = merged_counts
  )

  file_path <- file.path(outputDir, paste0("Table-", dv, "_Analysis_All.xlsx"))
  write_xlsx(output_list[[dv]], path = file_path)
}

## INTERACTION ANALYSIS for betaSupp_selfBase_pChange ##
## posthoc testing for LSM model
m_LSM <- lmer(betaSupp_selfBase_pChange ~ subgroup * cue + age + (1 | subject), df_LSM)
m_LSM_vcov_cr <- vcovCR(m_LSM, type = "CR2", cluster = model.frame(m_LSM)$subject)
m_LSM_result <- car::Anova(m_LSM, type = 3, vcov = m_LSM_vcov_cr)
emm_interaction1 <- emmeans(m_LSM, ~ subgroup * cue)
contrast_results1 <- club_contrast(m_LSM, m_LSM_vcov_cr, emm_interaction1, inter = TRUE)
contrast_df1 <- as.data.frame(contrast_results1) %>%
  mutate(
    p_val.fdr = p.adjust(p_val, method = "fdr"),
    signif.fdr = case_when(
      p_val.fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## posthoc testing for RSM model
m_RSM <- lmer(betaSupp_selfBase_pChange ~ subgroup * cue + age + (1 | subject), df_RSM)
m_RSM_vcov_cr <- vcovCR(m_RSM, type = "CR2", cluster = model.frame(m_RSM)$subject)
m_RSM_result <- car::Anova(m_RSM, type = 3, vcov = m_RSM_vcov_cr)
emm_interaction2 <- emmeans(m_RSM, ~ subgroup * cue)
contrast_results2 <- club_contrast(m_RSM, m_RSM_vcov_cr, emm_interaction2, inter = TRUE)
contrast_df2 <- as.data.frame(contrast_results2) %>%
  mutate(
    p_val.fdr = p.adjust(p_val, method = "fdr"),
    signif.fdr = case_when(
      p_val.fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

contrast_all <- rbind(
  contrast_df1 %>% mutate(model = "LSM"),
  contrast_df2 %>% mutate(model = "RSM")
)

contrast_all <- calculate_tstat(contrast_all, "Est", "SE") %>%
  mutate(across(where(is.character), ~ replace_na(., "")))

file_path <- file.path(outputDir, "Table-BetaSupp_Interaction_Contrasts.csv")
write.csv(contrast_all, file_path, row.names = FALSE)

################################################
## clean up table for for supplementary table ##
################################################
results_cue_pair_combined_clean <- output_list$betaSupp_selfBase_pChange$Cue_Pair
results_cue_pair_combined_clean <- results_cue_pair_combined_clean %>%
  select(variable, subgroup, contrast, Est, SE, df, CI_L, CI_U, p_val, p_val.fdr)

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
    Est = ifelse(flip, -Est, Est),
    CI_L = ifelse(flip, -orig_CI_U, orig_CI_L),
    CI_U = ifelse(flip, -orig_CI_L, orig_CI_U)
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
  CI95_Lower = CI_L,
  CI95_Upper = CI_U
)

results_cue_pair_combined_clean$Variable <- dplyr::recode(
  results_cue_pair_combined_clean$Variable,
  "LSM betaSupp_selfBase_pChange" = "LSM Contralateral Swing Beta (%)",
  "RSM betaSupp_selfBase_pChange" = "RSM Contralateral Swing Beta (%)"
)

results_cue_pair_combined_clean$Subgroup <- dplyr::recode(
  results_cue_pair_combined_clean$Subgroup,
  "hc" = "HC",
  "nofog" = "PD-NF",
  "fog" = "PD-F"
)

results_cue_pair_combined_clean <- calculate_tstat(results_cue_pair_combined_clean, "Est", "SE") %>%
  mutate(across(where(is.character), ~ replace_na(., ""))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

file_path <- file.path(outputDir, "Table-BetaSupp_Pairwise_Contrasts_ContralateralBeta.csv")
write.csv(results_cue_pair_combined_clean, file_path, row.names = FALSE)

warnings()
print('done.')
