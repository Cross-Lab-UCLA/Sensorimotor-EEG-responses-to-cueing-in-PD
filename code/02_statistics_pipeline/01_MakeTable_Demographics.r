## Make Demographics Table
# Summarizes demographic and clinical characteristics by subgroup (HC, PD-NF, PD-F).
# Performs ANOVA/Fisher's exact tests for group differences and saves results
#
# LM 120525
##

rm(list = ls())
library(tidyverse)
library(emmeans)
library(broom)
library(dplyr)
library(readxl)
library(gt)
library(openxlsx)
library(rmarkdown)

## set up paths
currentDir <- "E:/clab/DoD-Gait/results" # <-- Change this to your local path

dataDir <- file.path(currentDir, "processed_data")
outputDir <- file.path(currentDir, "output")

## load in demo
demo <- read_excel(file.path(dataDir, "demo_data_analysis.xlsx"))
demo <- demo %>%
  mutate(subgroup = fct_recode(subgroup,
                               "HC"     = "hc",
                               "PD-NF"  = "nofog",
                               "PD-F"    = "fog"))


demo <- demo %>%
  mutate(subgroup = factor(subgroup, levels = c("HC", "PD-NF", "PD-F")))
subgroup_counts <- demo %>%
  count(subgroup) %>%
  mutate(label = paste0(subgroup, " (n=", n, ")"))
levels(demo$subgroup) <- subgroup_counts$label

demo_sum <- demo %>%
  group_by(subgroup) %>%
  summarise(
    metronome_increased_sum = sum(metronome_increased, na.rm = TRUE),
    prop = mean(metronome_increased, na.rm = TRUE),
    n = n()
  )

# run ANOVA
vars <- c("Height", "Weight", "Age", "BENTON", "MOCA", "UPDRS3", "PAS","TrailBA")
#demo_anova <- vars %>%
#  set_names() %>%
#  map_df(~ {
#    formula <- as.formula(paste(.x, "~ subgroup"))
#    aov_result <- aov(formula, data = demo)
#    tidy(aov_result) 
#  }, .id = "variable")

results_list <- vars %>%
  set_names() %>%
  map(~ {
    formula <- as.formula(paste(.x, "~ subgroup"))
    model <- aov(formula, data = demo)
    anova_result <- tidy(model)
    
    # Check if the subgroup effect is significant
    p_val <- anova_result %>%
      filter(term == "subgroup") %>%
      pull(p.value)
    
    if (!is.null(p_val) && p_val < 0.05) {
      # Post hoc test
      emms <- emmeans(model, "subgroup")
      posthoc <- pairs(emms) %>%
        as.data.frame() %>%
        mutate(variable = .x)
    } else {
      posthoc <- NULL
    }
    
    list(anova = anova_result %>% mutate(variable = .x), posthoc = posthoc)
  })

demo_anova <- bind_rows(map(results_list, "anova"))
demo_posthoc <- bind_rows(map(results_list, "posthoc"))

demo_summary <- demo %>%
  group_by(subgroup) %>%
  summarise(
    across(
      .cols = all_of(vars),
      .fns = list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    -subgroup,
    names_to = c("variable", ".value"),
    names_sep = "_"
  ) %>%
  mutate(
    stat = sprintf("%.2f ± %.2f", mean, sd)
  ) %>%
  select(subgroup, variable, stat) %>%
  pivot_wider(names_from = subgroup, values_from = stat)


sex_summary <- demo %>%
  group_by(subgroup, Sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Sex, values_from = n, values_fill = 0) %>%
  mutate(stat = paste0("F:M: ", F, ":", M)) %>%
  select(subgroup, stat) %>%
  mutate(variable = "Sex") %>%
  pivot_wider(names_from = subgroup, values_from = stat)

affectedSide_summary <- demo %>%
  group_by(subgroup, Affected_Side) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Affected_Side, values_from = n, values_fill = 0) %>%
  mutate(stat = paste0("L:R: ", L, ":", R)) %>%
  select(subgroup, stat) %>%
  mutate(variable = "Affected_Side") %>%
  pivot_wider(names_from = subgroup, values_from = stat)

# run Fisher
#demo <- demo %>%
#  filter(subgroup != "hc")
#tbl <- table(demo[["Affected_Side"]], demo$subgroup)
#fisher.test(tbl)

vars_f <- c("Sex","Affected_Side")
demo_fisher <- vars_f %>%
  set_names() %>%
  map_df(~ {
    tbl <- table(demo[[.x]], demo$subgroup)
    print(.x)
    print(tbl)

    test <- fisher.test(tbl)
    tidy(test) %>%
      mutate(variable = .x)
  }, .id = NULL) %>%
  select(variable, p.value, method)

fog_summary <- demo %>%  
  group_by(subgroup) %>%
  summarise(
    across(
      .cols = c(FOG_Q),
      .fns = list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    stat = sprintf("%.2f ± %.2f", FOG_Q_mean, FOG_Q_sd),
    variable = "FOG_Q"
  ) %>%
  select(variable, subgroup, stat) %>%
  pivot_wider(names_from = subgroup, values_from = stat)

demo_summary <- bind_rows(sex_summary, affectedSide_summary, fog_summary, demo_summary)

anova_summary <- demo_anova %>%
  filter(term == "subgroup") %>%
  select(variable, p.value) %>%
  mutate(
    p.value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )

fisher_summary <- demo_fisher %>%
  select(variable, p.value) %>%
  mutate(
    p.value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )

bind_tests <- bind_rows(anova_summary, fisher_summary)

final_table <- demo_summary %>%
  left_join(bind_tests, by = "variable")

gt_table <- final_table %>%
  gt() %>%
  tab_header(title = "Demographics") %>%
  fmt_number(columns = where(is.numeric), decimals = 3) %>%
  cols_label(
    variable = " ",
    p.value = "p-value"
  )

# save as png
gtsave(gt_table, filename = "Table-demographics.png", path = outputDir)

# save as excel
wb <- createWorkbook()
addWorksheet(wb, "Demographics")
writeData(wb, "Demographics", gt_table)
addWorksheet(wb, "Post-Hoc")
writeData(wb, "Post-Hoc", demo_posthoc)
saveWorkbook(wb, file = file.path(outputDir, "Table-demographics.xlsx"), overwrite = TRUE)
demo_posthoc
