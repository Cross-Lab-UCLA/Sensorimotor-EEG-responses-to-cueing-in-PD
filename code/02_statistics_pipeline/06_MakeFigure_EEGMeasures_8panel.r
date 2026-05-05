## Make beta plots for manuscript
#
# LM 050426
##

# set up libraries
rm(list = ls())
library(lmerTest)
library(clubSandwich)
library(tidyverse)
library(readxl)
library(emmeans)
library(sjPlot)
library(ggeffects)
library(ggnewscale)
library(performance)
library(ggsignif)
library(patchwork)
library(scales)


# Set paths
currentDir <- "E:/clab/DoD-Gait/results" # <-- Change this to your local path

dataDir <- file.path(currentDir, "processed_data")
outputDir <- file.path(currentDir, "output")
funcDir <- file.path(currentDir,"a1_statistics_pipeline","func")

# load in functions
source(file.path(funcDir, "club_contrast.r"))
source(file.path(funcDir, "calculate_tstat.r"))

# load p values from stat tables
cue_pval_1 <- read.csv(file.path(
  outputDir,
  "Table-betaSupp_selfBase_pChange_4plot_Cue_Pairwise.csv"
))
cue_pval_affectedSide <- read.csv(file.path(
  outputDir,
  "Table-betaSupp_selfBase_pChange_4plot_Cue_Pairwise_affectedSide.csv"
))
cue_pval_all <- bind_rows(
  cue_pval_1,
  cue_pval_affectedSide
)

# set variables for plotting
dv <- "betaSupp_selfBase_pChange"
files <- c("RSM_erspData.csv", "LSM_erspData.csv", "moreAffected_erspData.csv", "lessAffected_erspData.csv")

#####################
### plot boxplots ###
#####################
plot_box <- function(filename, dv, cue_pval) {
  # Derive a readable ROI name (e.g., RSM, LSM, etc.)
  roi <- str_remove(filename, "_erspData.csv")
  title <- case_when(
    roi == "RSM" ~ "Right Sensorimotor",
    roi == "LSM" ~ "Left Sensorimotor",
    roi == "moreAffected" ~ "More Affected Hemisphere",
    roi == "lessAffected" ~ "Less Affected Hemisphere",
    TRUE ~ roi
  )
  df <- read_csv(file.path(dataDir, filename))

  df <- df %>%
    mutate(
      subgroup = factor(subgroup, levels = c("hc", "nofog", "fog")),
      cue = dplyr::recode(cue,
        "nocue" = "No Cue",
        "auditory" = "Auditory",
        "visual" = "Visual"
      ),
      cue = factor(cue, levels = c("Auditory", "No Cue", "Visual")),
      subgroup_jit = jitter(
        as.numeric(factor(subgroup)) + 0.75 / 3 *
          ifelse(cue == "Auditory", -1, ifelse(cue == "No Cue", 0, 1)),
        amount = 0.03
      )
    )

  sample_sizes <- df %>%
    group_by(subgroup) %>%
    summarise(n = n_distinct(subject), .groups = "drop") %>%
    mutate(label = paste0(
      case_when(
        subgroup == "hc" ~ "HC",
        subgroup == "nofog" ~ "PD-NF",
        subgroup == "fog" ~ "PD-F",
        TRUE ~ as.character(subgroup)
      ),
      "\n(n=", n, ")"
    )) %>%
    pull(label, subgroup)

  p <- ggplot(df, aes(x = subgroup, y = !!sym(dv), col = cue)) +
    geom_boxplot(alpha = 0.01, size = 1, outlier.shape = NA) +
    geom_point(aes(x = subgroup_jit, y = !!sym(dv)), alpha = 0.6) +
      scale_color_manual(values = c("#FF8C00", "#bd4646", "#0072b2")) +
    labs(title = title, x = NULL, y = "\n Contralateral Beta %") +
    theme_minimal(base_size = 18) +
    guides(color = guide_legend(title = "Cue Type")) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
    scale_x_discrete(labels = sample_sizes) +
    new_scale_color() +
    geom_line(
      aes(x = subgroup_jit, y = !!sym(dv), group = subject, color = cue),
      alpha = 0.4, size = 0.9, show.legend = FALSE
    ) +
      scale_color_manual(values = c("#FF8C00", "#0072b2", "#bd4646"))

  # add and plot significant bar
  ann_df <- cue_pval_all %>%
    filter(variable == paste0(roi, " ", dv), grepl("\\+", signif))

  if (nrow(ann_df) > 0) {
    ymax <- max(df[[dv]], na.rm = TRUE)
    ann_df <- ann_df %>%
      arrange(subgroup) %>%
      group_by(subgroup) %>%
      mutate(
        y = ymax + 1 + cumsum(rep(3, n())),
      ) %>%
      ungroup()

    ann_corrected <- ann_df %>% filter(!is.na(signif.fdr) & signif.fdr != "")
    ann_uncorrected <- ann_df %>% filter(is.na(signif.fdr) | signif.fdr == "")

    p <- p + geom_signif(
      inherit.aes = FALSE,
      data = ann_corrected,
      aes(
        xmin = min, xmax = max, annotations = signif.fdr, y_position = y, group = Est,
      ),
      manual = TRUE, tip_length = 0.02, textsize = 7.5, vjust = 0.5,
      color = "black"
    ) + geom_signif(
      inherit.aes = FALSE,
      data = ann_uncorrected,
      aes(
        xmin = min, xmax = max, annotations = signif, y_position = y, group = Est,
      ),
      manual = TRUE, tip_length = 0.02, textsize = 6, vjust = 0.1,
      color = "grey10"
    )
  }
  p
}

plots_box <- map(files, ~ plot_box(.x, dv = dv, cue_pval = cue_pval_all))
names(plots_box) <- str_remove(files, "_erspData.csv") # Name by ROI

######################
## plot interaction ##
######################
plot_int <- function(filename, dv, cue_pval) {
  roi <- str_remove(filename, "_erspData.csv")
  title <- case_when(
    roi == "RSM" ~ "Right Sensorimotor",
    roi == "LSM" ~ "Left Sensorimotor",
    roi == "moreAffected" ~ "More Affected Side",
    roi == "lessAffected" ~ "Less Affected Side",
    TRUE ~ roi
  )

  df <- read_csv(file.path(dataDir, filename)) %>% as.data.frame()
  df <- df %>%
    mutate(
      subgroup = factor(subgroup, levels = c("hc", "nofog", "fog")),
      cue = dplyr::recode(cue, "nocue" = "No Cue", "auditory" = "Auditory", "visual" = "Visual"),
      cue = factor(cue, levels = c("Auditory", "No Cue", "Visual")),
      subgroup_jit = jitter(
        as.numeric(subgroup) + 0.75 / 3 * ifelse(cue == "Auditory", -1, ifelse(cue == "No Cue", 0, 1)),
        amount = 0.03
      )
    )

  sample_sizes <- df %>%
    group_by(subgroup) %>%
    summarise(n = n_distinct(subject), .groups = "drop") %>%
    mutate(label = paste0(
      case_when(
        subgroup == "hc" ~ "HC",
        subgroup == "nofog" ~ "PD-NF",
        subgroup == "fog" ~ "PD-F",
        TRUE ~ as.character(subgroup)
      ),
      "\n(n=", n, ")"
    )) %>%
    pull(label, subgroup)

  form <- as.formula(paste0(dv, " ~ subgroup * cue + age + (1 | subject)"))
  m1 <- lmer(form, data = df)
  m1_vcov_cr <- vcovCR(m1, type = "CR2", cluster = model.frame(m1)$subject)

  emm_int <- emmeans(m1, ~ subgroup * cue)
  m1c_int <- club_contrast(m1, m1_vcov_cr, emm_int, TRUE) %>%
    mutate(
      p_val.fdr = p.adjust(p_val, "fdr"),
      signif.fdr = case_when(
        p_val.fdr < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  # Cue vs No Cue contrasts
  cue_vs_nocue <- emmeans(m1, ~ cue | subgroup) %>%
    club_contrast(m1, m1_vcov_cr, ., FALSE) %>%
    filter(str_detect(contrast, "No Cue")) %>%
    mutate(
      cue1 = str_trim(str_extract(contrast, "^[^-]+")),
      cue2 = str_trim(str_extract(contrast, "[^-]+$")),
      contrast = if_else(cue1 == "No Cue", paste(cue2, "-", cue1), contrast),
      across(c(Est, CI_L, CI_U), ~ if_else(cue1 == "No Cue", -.x, .x))
    )

  # Base plot
  p <- ggplot(cue_vs_nocue, aes(x = subgroup, y = Est, color = contrast, group = contrast)) +
    geom_point(position = position_dodge(0.4), size = 4) +
    geom_errorbar(aes(ymin = CI_L, ymax = CI_U),
      position = position_dodge(0.4), width = 0.25
    ) +
    scale_x_discrete(labels = sample_sizes) +
    labs(
      title = paste0(title, " Interaction"),
      x = NULL,
      y = paste0("Cue - No-Cue \n Contralateral Beta Δ% \n "),
      color = "Cue Contrast"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(colour = "black")
    ) +
      scale_color_manual(values = c("#FF8C00", "#0072b2"))

  # Significance bars
  max_y <- max(cue_vs_nocue$CI_U, na.rm = TRUE)
  has_hc <- any(str_detect(m1c_int$subgroup_pairwise, "hc"))
  group_map <- if (has_hc) {
    c(hc = 1, nofog = 2, fog = 3)
  } else {
    c(nofog = 1, fog = 2)
  }

  sig_contrasts1 <- m1c_int %>%
    filter(str_detect(cue_pairwise, "No Cue"), str_detect(signif, "\\+")) %>%
    mutate(
      left_grp = str_extract(subgroup_pairwise, "^[A-Za-z]+") %>% trimws(),
      right_grp = str_extract(subgroup_pairwise, "(?<=-)\\s*[A-Za-z]+$") %>% trimws(),
      xmin = group_map[left_grp],
      xmax = group_map[right_grp],
      shift = case_when(
        str_detect(cue_pairwise, regex("Auditory", ignore_case = TRUE)) ~ -0.1,
        str_detect(cue_pairwise, regex("Visual", ignore_case = TRUE)) ~ 0.1,
        TRUE ~ 0
      ),
      xmin = xmin + shift,
      xmax = xmax + shift,
      y = 5 + seq_along(xmin)
    )

  if (nrow(sig_contrasts1) > 0) {
    ann_corrected <- sig_contrasts1 %>% filter(!is.na(signif.fdr) & signif.fdr != "")
    ann_uncorrected <- sig_contrasts1 %>% filter(is.na(signif.fdr) | signif.fdr == "")

    p <- p + geom_signif(
      data = ann_corrected,
      aes(xmin = xmin, xmax = xmax, annotations = signif.fdr, y_position = y, group = Est),
      manual = TRUE, inherit.aes = FALSE, tip_length = 0.03, textsize = 8.5, vjust = 0.4,
      color = "black", textcolor = "black"
    ) + geom_signif(
      data = ann_uncorrected,
      aes(xmin = xmin, xmax = xmax, annotations = signif, y_position = y, group = Est),
      manual = TRUE, inherit.aes = FALSE, tip_length = 0.03, textsize = 6.5, vjust = 0.1,
      color = "grey10"
    )
  }
  p
}

plots_int <- map(files, ~ plot_int(.x, dv = dv, cue_pval = cue_pval_all))
names(plots_int) <- str_remove(files, "_erspData.csv") # Name by ROI

##################
### Save plots ###
##################

# set up boxplots
all_data <- map(plots_box, "data")
y_min <- min(map_dbl(all_data, ~ min(.x$betaSupp_selfBase_pChange, na.rm = TRUE)))
y_max <- max(map_dbl(all_data, ~ max(.x$betaSupp_selfBase_pChange, na.rm = TRUE)))

ylim_vals <- c(-50, 25)

p_lsm <- plots_box$LSM + coord_cartesian(xlim = c(0.5, 3.5), ylim = ylim_vals, expand = FALSE) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), ,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1)

p_rsm <- plots_box$RSM + coord_cartesian(xlim = c(0.5, 3.5), ylim = ylim_vals, expand = FALSE) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), ,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1)

p_more <- plots_box$moreAffected + coord_cartesian(xlim = c(0.5, 2.5), ylim = ylim_vals, expand = FALSE) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), ,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1)

p_less <- plots_box$lessAffected + coord_cartesian(xlim = c(0.5, 2.5), ylim = ylim_vals, expand = FALSE) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), ,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1)

# set up int plots
ylim_vals <- c(-12, 12)
p2_lsm <- plots_int$LSM + coord_cartesian(xlim = c(0.5, 3.5), ylim = ylim_vals, expand = FALSE) +
  scale_y_continuous(breaks = seq(ylim_vals[1], ylim_vals[2], length.out = 5)) +
  theme(
    plot.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), ,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1)

p2_rsm <- plots_int$RSM + coord_cartesian(xlim = c(0.5, 3.5), ylim = ylim_vals, expand = FALSE) +
  scale_y_continuous(breaks = seq(ylim_vals[1], ylim_vals[2], length.out = 5)) +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), ,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1)

p2_more <- plots_int$moreAffected + coord_cartesian(xlim = c(0.5, 2.5), ylim = ylim_vals, expand = FALSE) +
  scale_y_continuous(breaks = seq(ylim_vals[1], ylim_vals[2], length.out = 5)) +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), ,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1)

p2_less <- plots_int$lessAffected + coord_cartesian(xlim = c(0.5, 2.5), ylim = ylim_vals, expand = FALSE) +
  scale_y_continuous(breaks = seq(ylim_vals[1], ylim_vals[2], length.out = 5)) +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), ,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1)

# combine plots
c_plot <- (p_lsm | p_rsm | p_more | p_less) /
  (p2_lsm | p2_rsm | p2_more | p2_less) +
  plot_layout(guides = "collect") &
  theme(
    text = element_text(size = 18),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.text = element_text(size = 18)
  )

# Save plot
saveFile <- file.path(outputDir, "Figure-BetaSupp_Contrast_8Panel.jpg")
ggsave(saveFile, plot = c_plot, width = 18, height = 9, dpi = 600)

print("done.")
