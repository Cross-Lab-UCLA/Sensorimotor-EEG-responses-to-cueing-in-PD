## Make stat tables for the more and less side models
# dv ~ cue x group + age + pas + updrs3 + (1 | subject)
# Require data from A10_Extract2CSV_affectedSide_mean
#
# LM 050426
##

# set up libraries
rm(list = ls())
library(lmerTest)
library(car)
library(tidyverse)
library(readxl)
library(writexl)
library(emmeans)
library(ggeffects)
library(performance)
library(ggsignif)
library(clubSandwich)

# set up paths
currentDir <- "E:/clab/DoD-Gait/results" # <-- Change this to your local path

dataDir <- file.path(currentDir, "processed_data")
outputDir <- file.path(currentDir, "output")
funcDir <- file.path(currentDir,"a1_statistics_pipeline","func")

# load in functions
source(file.path(funcDir, "club_contrast.r"))
source(file.path(funcDir, "calculate_tstat.r"))

# dependent variables
dvs <- c("betaSupp_selfBase_pChange")

# filenames
MORE_filename <- "moreAffected_erspData.csv"
LESS_filename <- "lessAffected_erspData.csv"

drop_subjects <- c("")

clean_df <- function(filename, dataDir, drop_subjects) {
  df <- read_csv(file.path(dataDir, filename), show_col_types = FALSE)

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

df_MORE <- clean_df(MORE_filename, dataDir, drop_subjects)
df_LESS <- clean_df(LESS_filename, dataDir, drop_subjects)
output_list <- list()

# run loop
for (dv in dvs) {
  message("Running analysis for: ", dv)

  # set formula
  formula1 <- as.formula(paste(dv, "~ cue * subgroup + age + pas + updrs3 + (1 | subject)"))
  options(contrasts = c("contr.sum", "contr.sum"))

  # MORE affected side MODEL
  subgroup_counts <- df_MORE %>%
    distinct(subject, subgroup) %>%
    count(subgroup)
  for (i in seq_len(nrow(subgroup_counts))) {
    message(
      "More Affected Subgroup '", subgroup_counts$subgroup[i],
      "', n = ", subgroup_counts$n[i], " subjects."
    )
  }

  MORE <- lmer(formula1, data = df_MORE)
  MORE_vcov_cr <- vcovCR(MORE, type = "CR2", cluster = model.frame(MORE)$subject)
  MORE_a <- car::Anova(MORE, type = 3, vcov = MORE_vcov_cr)

  MORE_anova_df <- MORE_a %>%
    as.data.frame() %>%
    mutate(model = paste0("MORE ", dv)) %>%
    tibble::rownames_to_column(var = "term")

  emm_cueCompByGroup <- emmeans(MORE, ~ cue | subgroup)
  MORE_cue_pairs_by_subgroup <- club_contrast(MORE, MORE_vcov_cr, emm_cueCompByGroup) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  MORE_emm_cue <- emmeans(MORE, ~cue)
  MORE_cue_main <- club_contrast(MORE, MORE_vcov_cr, MORE_emm_cue) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  MORE_emm_group <- emmeans(MORE, ~subgroup)
  MORE_group_main <- club_contrast(MORE, MORE_vcov_cr, MORE_emm_group) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  # LESS MODEL
  subgroup_counts <- df_LESS %>%
    distinct(subject, subgroup) %>%
    count(subgroup)
  for (i in seq_len(nrow(subgroup_counts))) {
    message(
      "LESS Subgroup '", subgroup_counts$subgroup[i],
      "', n = ", subgroup_counts$n[i], " subjects."
    )
  }

  LESS <- lmer(formula1, data = df_LESS)
  LESS_vcov_cr <- vcovCR(LESS, type = "CR2", cluster = model.frame(LESS)$subject)
  LESS_a <- car::Anova(LESS, type = 3, vcov = LESS_vcov_cr)

  LESS_anova_df <- LESS_a %>%
    as.data.frame() %>%
    mutate(model = paste0("LESS ", dv)) %>%
    tibble::rownames_to_column(var = "term")

  emm_cueCompByGroup <- emmeans(LESS, ~ cue | subgroup)
  LESS_cue_pairs_by_subgroup <- club_contrast(LESS, LESS_vcov_cr, emm_cueCompByGroup) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  LESS_emm_cue <- emmeans(LESS, ~cue)
  LESS_cue_main <- club_contrast(LESS, LESS_vcov_cr, LESS_emm_cue) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  LESS_emm_group <- emmeans(LESS, ~subgroup)
  LESS_group_main <- club_contrast(LESS, LESS_vcov_cr, LESS_emm_group) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, method = "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variable = dv
    ) %>%
    select(variable, everything())

  # combine anova tables from LESS andMORE
  anova_combined <- bind_rows(LESS_anova_df, MORE_anova_df)
  anova_combined <- anova_combined %>%
    select(model, everything()) %>%
    mutate(
      signif = case_when(
        `Pr(>Chisq)` < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  ## planned contrasts ##
  cue_pairwise_combined <- bind_rows(
    LESS_cue_pairs_by_subgroup %>% mutate(variable = paste("lessAffected", dv)),
    MORE_cue_pairs_by_subgroup %>% mutate(variable = paste("moreAffected", dv))
  )

  cue_pairwise_combined <- cue_pairwise_combined %>%
    extract(
      contrast,
      into = c("comp1", "comp2"),
      regex = "([a-zA-Z ]+) - (.+)",
      remove = FALSE
    )

  # create table for contrasts that are signficant for plotting
  results_p <- cue_pairwise_combined %>%
    filter(p_val < 0.05) %>%
    mutate(
      group_val = case_when(
        subgroup == "nofog" ~ 1,
        subgroup == "fog" ~ 2,
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
    select(variable, subgroup, comp1, comp2, Est, SE, CI_L, CI_U, p_val, signif, p_val.fdr, signif.fdr, min, max)

  file_path <- file.path(outputDir, paste0("Table-", dv, "_4plot_Cue_Pairwise_affectedSide.csv"))
  write.csv(results_p, file_path, row.names = FALSE)

  ## main effects
  cue_main_combined <- bind_rows(
    LESS_cue_main %>% mutate(variable = paste("LESS", dv)),
    MORE_cue_main %>% mutate(variable = paste("MORE", dv))
  )

  group_main_combined <- bind_rows(
    LESS_group_main %>% mutate(variable = paste("LESS", dv)),
    MORE_group_main %>% mutate(variable = paste("MORE", dv))
  )

  # combine emmean cue
  cue_combined <- bind_rows(
    as.data.frame(MORE_emm_cue) %>% mutate(variable = paste("MORE", dv)),
    as.data.frame(LESS_emm_cue) %>% mutate(variable = paste("LESS", dv))
  ) %>%
    select(variable, everything())

  group_combined <- bind_rows(
    as.data.frame(MORE_emm_group) %>% mutate(variable = paste("MORE", dv)),
    as.data.frame(LESS_emm_group) %>% mutate(variable = paste("LESS", dv))
  ) %>%
    select(variable, everything())

  ### add t-stat
  cue_main_combined <- calculate_tstat(cue_main_combined, "Est", "SE") %>%
    mutate(across(where(is.character), ~ replace_na(., "")))

  group_main_combined <- calculate_tstat(group_main_combined, "Est", "SE") %>%
    mutate(across(where(is.character), ~ replace_na(., "")))

  cue_pairwise_combined <- calculate_tstat(cue_pairwise_combined, "Est", "SE") %>%
    mutate(across(where(is.character), ~ replace_na(., "")))

  ## save to excel
  output_list[[dv]] <- list(
    ANOVA = anova_combined,
    Cue_Est = cue_combined,
    Group_Est = group_combined,
    Cue_Main = cue_main_combined,
    Cue_Pair = cue_pairwise_combined,
    Group_Main = group_main_combined
  )

  file_path <- file.path(outputDir, paste0("Table-", dv, "_Analysis_All_affectedSide.xlsx"))
  write_xlsx(output_list[[dv]], path = file_path)
}
warnings()


## INTERACTION ANALYSIS for betaSupp_selfBase_pChange##
## posthoc testing for MORE model
m_MORE <- lmer(betaSupp_selfBase_pChange ~ subgroup * cue + age + pas + updrs3 + (1 | subject), df_MORE)
m_MORE_vcov_cr <- vcovCR(m_MORE, type = "CR2", cluster = model.frame(m_MORE)$subject)
m_MORE_result <- car::Anova(m_MORE, type = 3, vcov = m_MORE_vcov_cr)
emm_interaction1 <- emmeans(m_MORE, ~ subgroup * cue)
contrast_results1 <- club_contrast(m_MORE, m_MORE_vcov_cr, emm_interaction1, inter = TRUE)
contrast_df1 <- as.data.frame(contrast_results1) %>%
  mutate(
    p_val.fdr = p.adjust(p_val, method = "fdr"),
    signif.fdr = case_when(
      p_val.fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

## posthoc testing for LESS model
m_LESS <- lmer(betaSupp_selfBase_pChange ~ subgroup * cue + age + pas + updrs3 + (1 | subject), df_LESS)
m_LESS_vcov_cr <- vcovCR(m_LESS, type = "CR2", cluster = model.frame(m_LESS)$subject)
m_LESS_result <- car::Anova(m_LESS, type = 3, vcov = m_LESS_vcov_cr)
emm_interaction2 <- emmeans(m_LESS, ~ subgroup * cue)
contrast_results2 <- club_contrast(m_LESS, m_LESS_vcov_cr, emm_interaction2, inter = TRUE)
contrast_df2 <- as.data.frame(contrast_results2) %>%
  mutate(
    p_val.fdr = p.adjust(p_val, method = "fdr"),
    signif.fdr = case_when(
      p_val.fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

contrast_all <- rbind(
  contrast_df1 %>% mutate(model = "MORE"),
  contrast_df2 %>% mutate(model = "LESS")
)

contrast_all <- calculate_tstat(contrast_all,"Est", "SE") %>%
  mutate(across(where(is.character), ~ replace_na(., "")))

file_path <- file.path(outputDir, "Table-BetaSupp_Interaction_Contrasts_affectedSide.csv")
write.csv(contrast_all, file_path, row.names = FALSE)
warnings()

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
  "lessAffected betaSupp_selfBase_pChange" = "Less Affected Side Contralateral Beta (%)",
  "moreAffected betaSupp_selfBase_pChange" = "More Affected Side Contralateral Beta (%)"
)

results_cue_pair_combined_clean$Subgroup <- dplyr::recode(
  results_cue_pair_combined_clean$Subgroup,
  "hc" = "HC",
  "nofog" = "PD-NF",
  "fog" = "PD-F"
)

results_cue_pair_combined_clean <- calculate_tstat(results_cue_pair_combined_clean,"Est", "SE") %>%
  mutate(across(where(is.character), ~ replace_na(., ""))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

file_path <- file.path(outputDir, "Table-BetaSupp_Pairwise_Contrasts_ContralateralBeta_MoreAffected.csv")
write.csv(results_cue_pair_combined_clean, file_path, row.names = FALSE)

warnings()
print('done.')
