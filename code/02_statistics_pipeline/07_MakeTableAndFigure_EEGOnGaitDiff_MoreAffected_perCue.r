## Make stat tables for 2-way interaction regression model
# gait measure change ~ group * EEG measure change (per Cue)
#
# LM 050426
##

rm(list = ls())
library(tidyverse)
library(car)
library(emmeans)
library(glue)
library(patchwork)
library(writexl)
library(lmtest)
library(sandwich)
library(interactions)
library(jtools)
library(ggeffects)
library(ggplot2)

# set up paths
currentDir <- "E:/clab/DoD-Gait/results" # <-- Change this to your local path

dataDir <- file.path(currentDir, "processed_data")
outputDir <- file.path(currentDir, "output")
funcDir <- file.path(currentDir,"a1_statistics_pipeline","func")

# variables of interest
iv_list <- c(
  "betaSupp_selfBase_pChange_diffFromNoCue"
)

cols_to_exclude <- c("subject", "group", "subgroup", "cue", "age", "weight", "height", "moca", "updrs3", "pas")

#outcome_templates <- c(
#  Stride_mean = "mean_StrideLength_{limb}_diffFromNoCue",
#  Speed_mean = "mean_Speed_{limb}_diffFromNoCue",
#  Cadence_mean = "cv_Cadence_{limb}_diffFromNoCue",
#  Stride_cv = "cv_StrideLength_{limb}_diffFromNoCue",
#  Speed_cv = "cv_Speed_{limb}_diffFromNoCue",
#  Cadence_cv = "cv_Cadence_{limb}_diffFromNoCue"
#)

outcome_templates <- c(
  Stride_mean = "mean_StrideLength_{limb}_diffFromNoCue",
  Speed_mean = "mean_Speed_{limb}_diffFromNoCue",
  Cadence_cv = "cv_Cadence_{limb}_diffFromNoCue"
)

## set up model
# dv ~ cue * subgroup * iv + age + (1 | subject)
##
build_formulas <- function(iv, limb) {
  map(outcome_templates, ~ {
    outcome <- str_replace(.x, "\\{limb\\}", limb)
    #as.formula(glue("{outcome} ~ subgroup * {iv} + age + updrs3 + pas + (1 | subject)"))
    as.formula(glue("{outcome} ~ subgroup * {iv} + age + updrs3 + pas"))
  })
}

# load in data and perpare data for analysis, including
# recaling, substacing by baseline, factoring the levels
preprocess_data <- function(df, cols_to_exclude, iv, limb, baseline_cue = "nocue") {
  # check if the which side being processed
  if (grepl("affected", limb, ignore.case = TRUE)) {
    affectedSide <- TRUE
    message(glue("Processing data for {limb}"))
  } else {
    affectedSide <- FALSE
  }
  message(glue("Processing data for {limb}"))

  # check stride and speed variables for rescaling to cm
  matched_cols <- df %>%
    select(matches("(?i)(mean_stride|mean_speed)")) %>%
    names()
  message("Scaling (to cm) variables matched for stride/speed: ", paste(matched_cols, collapse = ", "))
  df <- df %>%
    mutate(across(matches("(?i)(mean_stride|mean_speed)"), ~ .x * 100)) %>%
    filter(subgroup != "hc")

  # check current levels of subgroup BEFORE recoding
  current_subgroups <- unique(as.character(df$subgroup))
  expected_levels <- c("fog", "nofog")

  missing_levels <- setdiff(expected_levels, current_subgroups)
  if (length(missing_levels) > 0) {
    warning(glue("Missing expected subgroup levels in data before recoding: {paste(missing_levels, collapse = ', ')}"))
    warning(glue("Current subgroup levels are: {paste(current_subgroups, collapse = ', ')}"))
  } else {
    message(glue("All expected subgroup levels present: {paste(expected_levels, collapse = ', ')}"))
  }

  ## subject by no cue condition
  cols_to_subtract <- df %>%
    select(-all_of(cols_to_exclude)) %>%
    select(where(is.numeric)) %>%
    names()

  # warnings if there is no nocue
  df <- df %>%
    group_by(subject, subgroup) %>%
    mutate(across(
      all_of(cols_to_subtract),
      ~ {
        if (!any(cue == baseline_cue, na.rm = TRUE)) {
          warning(glue("No '{baseline_cue}' condition for subject={unique(subject)} "))
          return(NA_real_)
        }
        baseline_vals <- .x[cue == baseline_cue]
        if (length(na.omit(baseline_vals)) > 1) {
          warning(glue("Multiple '{baseline_cue}' rows for subject={unique(subject)} — using first."))
        }
        .x - first(na.omit(baseline_vals))
      },
      .names = "{.col}_diffFromNoCue"
    )) %>%
    ungroup()

  # normalize stride length if column exists
  stride_col <- glue("mean_StrideLength_{limb}_diffFromNoCue")
  if (stride_col %in% names(df)) {
    df <- df %>%
      mutate(!!stride_col := !!sym(stride_col) / height)
    message(glue("Normalized {stride_col} by height."))
  } else {
    message(glue("Column {stride_col} not found — skipping normalization."))
  }

  df <- df %>%
    filter(cue != baseline_cue) %>%
    select(all_of(cols_to_exclude), ends_with("_diffFromNoCue")) %>%
    mutate(
      cue = fct_relevel(cue, "auditory"),
      subgroup = fct_recode(subgroup, "PD-F" = "fog", "PD-NF" = "nofog"),
      subgroup = fct_relevel(subgroup, "PD-NF", "PD-F")
    )

  # Set contrasts for modeling
  contrasts(df$cue) <- contr.sum
  contrasts(df$subgroup) <- contr.sum

  return(df)
}


# model fitting function
fit_model <- function(formula, data, iv, title) {
  
  model <- lm(formula, data = data)
  vcov_model <- vcovHC(model, type = "HC3")
  model_a <- car::Anova(model, type = 3, vcov = vcov_model)
  anova_tbl <- model_a %>%
    as.data.frame() %>%
    rownames_to_column("Predictors")
  simple_slopes <- as.data.frame(
    summary(
      emtrends(model, ~ subgroup, var = iv, vcov = vcov_model),
      infer = c(TRUE, TRUE) # <- must go here
    )
  )
  names(simple_slopes)[grepl("\\.trend$", names(simple_slopes))] <- "estimate"


  # Contrast of slopes
  contrast_slopes <- as.data.frame(
    summary(
      pairs(
        emtrends(model, ~subgroup, var = iv, vcov = vcov_model)
      ),
      infer = c(TRUE, TRUE) # <- and here
    )
  )

  simple_slopes$type <- "simple_slope"
  contrast_slopes$type <- "contrast"

  # Combine vertically
  all_slopes <- bind_rows(simple_slopes, contrast_slopes)
  all_slopes <- all_slopes %>%
    select(type, contrast, everything())

  list(model = model, anova = anova_tbl, slopes = all_slopes)
}

# Combine model summaries for output
combine_anova <- function(results, side_label) {
  imap_dfr(results, ~ .x$anova %>% mutate(Outcome = .y, Side = side_label, .before = 1))
}

combine_slopes <- function(results, side_label) {
  imap_dfr(results, ~ .x$all_slopes %>% mutate(Outcome = .y, Side = side_label, .before = 1))
}
cue_types <- c("visual", "auditory")

for (cue_type in cue_types) {
  message("Running cue type: ", cue_type)

  for (iv in iv_list) {

    message("\n\nRunning Analysis for: ", iv)
    iv_tag <- str_extract(iv, "^[a-zA-Z]+")

    df_MORE <- read_csv(file.path(dataDir, "moreAffected_erspData.csv"), show_col_types = FALSE) %>%
      preprocess_data(cols_to_exclude, iv, limb = "moreAffected") %>%
      filter(cue == cue_type)

    df_LESS <- read_csv(file.path(dataDir, "lessAffected_erspData.csv"), show_col_types = FALSE) %>%
      preprocess_data(cols_to_exclude, iv, limb = "lessAffected") %>%
      filter(cue == cue_type)

    formulas_MORE <- build_formulas(iv, "moreAffected")
    formulas_LESS <- build_formulas(iv, "lessAffected")

    message("Running MORE models...")
    MORE_results <- imap(formulas_MORE,
      ~ fit_model(.x, df_MORE, iv,
                  glue("MORE ({cue_type}) on {.y} Model"))
    )
MORE_results <- lapply(MORE_results, function(x) {
  x$slopes <- x$slopes %>% mutate(
    cue = cue_type,
    model_type = "MORE"
  ) %>%
    select(model_type, cue, everything())
    x
})

    message("Running LESS models...")
    LESS_results <- imap(formulas_LESS,
      ~ fit_model(.x, df_LESS, iv,
                  glue("LESS ({cue_type}) on {.y} Model"))
    )

LESS_results <- lapply(LESS_results, function(x) {
  x$slopes <- x$slopes %>%
    mutate(
      cue = cue_type,
      model_type = "LESS"
    ) %>%
    select(model_type, cue, everything())
    x
})


    all_anova <- bind_rows(
      combine_anova(MORE_results, "MORE"),
      combine_anova(LESS_results, "LESS")
    ) %>%
      mutate(signif = case_when(
        `Pr(>F)` < 0.001 ~ "***",
        `Pr(>F)` < 0.01 ~ "**",
        `Pr(>F)` < 0.05 ~ "*",
        TRUE ~ ""
      ))

    all_slopes <- bind_rows(
        MORE_results$Stride_mean$slopes %>% mutate(dv = "Stride"),
        LESS_results$Stride_mean$slopes %>% mutate(dv = "Stride"),
        MORE_results$Speed_mean$slopes %>% mutate(dv = "Speed"),
        LESS_results$Speed_mean$slopes %>% mutate(dv = "Speed"),
        MORE_results$Cadence_cv$slopes %>% mutate(dv = "Cadence"),
        LESS_results$Cadence_cv$slopes %>% mutate(dv = "Cadence")
    )

    ## fdr correction for slopes per dv and model family
    all_slopes <- all_slopes %>%
        group_by(dv, model_type) %>%
        mutate(
            p.fdr = p.adjust(p.value, method = "fdr")
        ) %>%
        ungroup()


    # SAVE FILES
    out_xlsx <- file.path(
      outputDir,
      glue("Table-BrainBehaviorReg_ANOVA_{iv_tag}_{cue_type}_clubSandwich.xlsx")
    )
    write_xlsx(
      list(ANOVA = all_anova, Slopes = all_slopes),
      path = out_xlsx
    )

    # Save significant-only
    sig_csv <- file.path(
      outputDir,
      glue("Table-BrainBehaviorReg_ANOVA_{iv_tag}_{cue_type}_clubSandwich_signifOnly.csv")
    )
    write_csv(all_anova %>% filter(signif != ""), sig_csv)

  } # end inner IV loop

} # end cue loop


##
## Plotting significant slopes for visual cue beta suppression on more affected speed change
##

df_MORE <- read_csv(file.path(dataDir, "moreAffected_erspData.csv"), show_col_types = FALSE) %>%
      preprocess_data(cols_to_exclude, iv, limb = "moreAffected") %>%
      filter(cue == "visual")

iv <- "betaSupp_selfBase_pChange_diffFromNoCue"
formula1 <- as.formula(
  paste("mean_Speed_moreAffected_diffFromNoCue ~  subgroup *", iv,
        "+ age + updrs3 + pas")
)

m1 <- lm(formula1, df_MORE)
#check_model(m1)
vcov_model1 <- vcovHC(m1, type = "HC3")
model_a <- car::Anova(m1, type = 3, vcov = vcov_model1)
pred_df <- ggpredict(
  m1,
  terms = c("betaSupp_selfBase_pChange_diffFromNoCue [-40:40 by=1]", "subgroup"),
  vcov = vcov_model1
)

p1 <- ggplot(pred_df, aes(
  x = x, y = predicted,
  linetype = group, color = group, fill = group
)) +
  geom_line(size = 2.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("#7300B3", "#0073B3")) +
  scale_fill_manual(values = c("#7300B3", "#0073B3")) +
  theme_minimal()
p1

p1 <- p1 +
  geom_point(
    data = df_MORE,
    inherit.aes = FALSE,
    aes(
      x = betaSupp_selfBase_pChange_diffFromNoCue,
      y = mean_Speed_moreAffected_diffFromNoCue,
      shape = subgroup,
      color = subgroup
    ),
    size = 6,
    alpha = 0.9
  ) +
  scale_color_manual(values = c("#7300B3", "#0073B3")) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) + 
  coord_cartesian(xlim = c(-25, 15), ylim = c(-30, 30)) +
  scale_x_continuous(breaks = seq(-25, 15, by = 10)) +
  scale_y_continuous(breaks = seq(-30, 30, by = 10)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1) + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1)

p1

counts <- df_MORE %>% 
  group_by(subgroup) %>% 
  summarise(n = n())

label_names <- setNames(
  paste0(counts$subgroup, " (n=", counts$n, ")"), 
  counts$subgroup
)

p1_save <- p1 + theme(
  text = element_text(size = 18, color = "black"),
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.justification = "center",
  axis.title = element_text(size = 28, color = "black", face = "plain"),
  axis.text = element_text(size = 22, color = "black"),
  legend.text = element_text(size = 22, color = "black"),
) + 
  labs(
    x = "Visual Cue - No-Cue Contralateral Beta (Δ %)",
    y = "Visual Cue - No-Cue Gait Speed (Δ cm/s)",
    color = "PD Group",
    shape = "PD Group",
    linetype = "PD Group"
  ) +
  scale_color_manual(name = "PD Group", values = c("#7300B3", "#0073B3"), labels = label_names) +
  scale_shape_manual(name = "PD Group", values = c(16, 17), labels = label_names) + 
  scale_fill_manual(name = "PD Group", values = c("#7300B3", "#0073B3"), labels = label_names) +
  scale_linetype_manual(name = "PD Group", values = c("dashed", "solid"), labels = label_names) +
  guides(fill = "none", linetype = "none")

p1_save

saveFile <- file.path(outputDir, "Figure-BetaSupp_Regression_1Panel.jpg")
ggsave(saveFile, plot = p1_save, width = 11, height = 9, dpi = 600)



##
## Plot 1 + Less Affected Cadence CV vs Auditory Cue Beta Suppression
##
df_LESS <- read_csv(file.path(dataDir, "lessAffected_erspData.csv"), show_col_types = FALSE) %>%
      preprocess_data(cols_to_exclude, iv, limb = "lessAffected") %>%
      filter(cue == "auditory")

formula2 <- as.formula(
  paste("cv_Cadence_lessAffected_diffFromNoCue ~  subgroup *", iv,
        "+ age + updrs3 + pas")
)

m2 <- lm(formula2, df_LESS)
#check_model(m2)
vcov_model2 <- vcovHC(m2, type = "HC3")
model_a <- car::Anova(m2, type = 3, vcov = vcov_model2)

coef_est <- coeftest(m2, vcov. = vcov_model2)["betaSupp_selfBase_pChange_diffFromNoCue", ]
coef_est_val <- coef_est["Estimate"]
se <- coef_est["Std. Error"]
ci <- coef_est_val + c(-1, 1) * qt(0.975, df.residual(m2)) * se
df_val <- df.residual(m2)

results_table <- data.frame(
  Term = "betaSupp_selfBase_pChange_diffFromNoCue",
  Estimate = coef_est_val,
  SE = se,
  Df = df_val,
  Lower.CL = ci[1],
  Upper.CL = ci[2],
  t.value = coef_est["t value"],
  p.value = coef_est["Pr(>|t|)"]
)

  save_csv <- file.path(
      outputDir,
      glue("Table-BrainBehaviorReg_MainCueEffect_slope.csv")
    )
  write_csv(results_table, save_csv)

p2 <- effect_plot(
  model = m2,
  pred = betaSupp_selfBase_pChange_diffFromNoCue,
  vcov = vcov_model2,
  data = df_LESS,
  interval = TRUE,
  plot.points = FALSE,
  point.size = 5,
  x.label = "Auditory Cue - No-Cue\n Contralateral Beta  (Δ %)",
  y.label = "Auditory Cue - No-Cue\n Cadence Variability (Δ %)",
)
p2 <- p2 +
  geom_point(
    data = df_LESS,
    aes(
      x = betaSupp_selfBase_pChange_diffFromNoCue,
      y = cv_Cadence_lessAffected_diffFromNoCue,
      shape = subgroup,
      color = subgroup
    ),
    size = 5,
    alpha = 0.9
  ) +
  scale_color_manual(values = c("#f7796d", "#9c0e02")) +
  labs(title = "Auditory Cue-Driven \nContralateral Beta and Cadence Variability Change \non Less-Affected Side") +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) + 
  coord_cartesian(xlim = c(-14.5, 7.35), ylim = c(-0.1, 0.1), expand = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1)
p2

c_plot <- p1 + p2
c_plot <- c_plot + plot_layout(guides = "collect") &
    theme(
        text = element_text(size = 18, color = "black"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        axis.title = element_text(size = 19, color = "black", face = "plain"),
        axis.text = element_text(size = 19, color = "black"),
        legend.text = element_text(size = 19, color = "black"),
        plot.title = element_text(size = 21, color = "black")
    )
c_plot

# Save plot
saveFile <- file.path(outputDir, "Figure-BetaSupp_Regression_2Panel.jpg")
ggsave(saveFile, plot = c_plot, width = 18, height = 9, dpi = 600)

