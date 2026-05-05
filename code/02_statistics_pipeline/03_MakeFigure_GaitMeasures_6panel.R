## Make kinematic plots for manuscript
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
source(file.path(currentDir, "club_contrast.r"))

# load p values
cue_pval <- read.csv(file.path(
    outputDir,
    "Table-Kinematics_4plot_Cue_Pairwise.csv"
))

# SET variables
dv_list <- c("speed_mean", "stride_mean", "cadence_cv")

df <- read_excel(file.path(dataDir, "steps_lmer_kinematic_run1.xlsx"))
df$group <- factor(df$group)
df$subgroup <- factor(df$subgroup, levels = c("hc", "nofog", "fog"))
df$subgroup <- relevel(df$subgroup, "hc")
df$group <- relevel(df$group, "hc")
df$cue <- gsub("nocue", "No Cue", df$cue) # change label
df$cue <- gsub("auditory", "Auditory", df$cue)
df$cue <- gsub("visual", "Visual", df$cue)
df$cue <- factor(df$cue, levels = unique(df$cue))
df$cue <- relevel(df$cue, "Auditory")
df <- df %>% # rescale stride and speed to cm
    mutate(
        stride_mean = stride_mean * 100 / height,
        speed_mean = speed_mean * 100
    )
df$subgroup_jit <- jitter(as.numeric(factor(df$subgroup)) + .9 / 3 *
    ifelse(df$cue == "Auditory", -1, ifelse(df$cue == "No Cue", 0, 1)), amount = .01)

## MUST SET CONSTRASTS for ANOVA
contrasts(df$cue) <- contr.sum
contrasts(df$group) <- contr.sum
contrasts(df$subgroup) <- contr.sum

# Calculate sample sizes for each subgroup
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

#####################
### plot boxplots ###
#####################
plot_box <- function(df, dv, cue_pval) {
    title <- case_when(
        dv == "speed_mean" ~ "Gait Speed",
        dv == "stride_mean" ~ "Stride Length",
        dv == "cadence_cv" ~ "Cadence Variability",
        TRUE ~ dv
    )

    ylabel <- case_when(
        dv == "speed_mean" ~ "\ncm/s",
        dv == "stride_mean" ~ "cm/height",
        dv == "cadence_cv" ~ "%",
        TRUE ~ dv
    )

    p <- ggplot(df, aes(x = subgroup, y = !!sym(dv), col = cue)) +
        geom_boxplot(alpha = 0.02, size = 1.2, outlier.shape = NA) +
        geom_point(aes(x = subgroup_jit, y = !!sym(dv)), alpha = 0.6) +
        scale_color_manual(values = c("#FF8C00", "#bd4646", "#0072b2")) +
        labs(title = title, x = NULL, y = ylabel) +
        theme_minimal(base_size = 18) +
        guides(color = guide_legend(title = "Cue Type")) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
        scale_x_discrete(labels = sample_sizes) +
        new_scale_color() +
        geom_line(
            aes(x = subgroup_jit, y = !!sym(dv), group = subject, color = cue),
            alpha = 0.4, linewidth = 0.9, show.legend = FALSE
        ) +
        scale_color_manual(values = c("#FF8C00", "#0072b2", "#bd4646"))

    # add and plot significant bar
    ann_df <- cue_pval %>%
        filter(variable == dv, grepl("\\+", signif))
    if (nrow(ann_df) > 0) {
        ymax <- max(df[[dv]], na.rm = TRUE)
        ann_df <- ann_df %>%
            arrange(subgroup) %>%
            group_by(subgroup) %>%
            mutate(y = (ymax * 0.99) + cumsum(rep(0.055 * ymax, n()))) %>%
            ungroup() %>%
            mutate(row_id = row_number())

        ann_corrected <- ann_df %>% filter(!is.na(signif.fdr) & signif.fdr != "")
        ann_uncorrected <- ann_df %>% filter(is.na(signif.fdr) | signif.fdr == "")

        p <- p + geom_signif(
            inherit.aes = FALSE,
            data = ann_corrected,
            aes(
                xmin = min, xmax = max, annotations = signif.fdr, y_position = y, group = Est,
            ),
            manual = TRUE, tip_length = 0.02, textsize = 7, vjust = 0.55,
            color = "black"
        ) + geom_signif(
            inherit.aes = FALSE,
            data = ann_uncorrected,
            aes(
                xmin = min, xmax = max, annotations = signif, y_position = y, group = Est,
            ),
            manual = TRUE, tip_length = 0.02, textsize = 6, vjust = 0.15,
            color = "grey10"
        )
    }
    p
}

plots_box <- map(dv_list, ~ plot_box(df, dv = .x, cue_pval = cue_pval))
names(plots_box) <- dv_list

######################
## plot interaction ##
######################
plot_int <- function(df, dv) {
    title <- case_when(
        dv == "speed_mean" ~ "Gait Speed",
        dv == "stride_mean" ~ "Stride Length",
        dv == "cadence_cv" ~ "Cadence Variability",
        TRUE ~ dv
    )

    ylabel <- case_when(
        dv == "speed_mean" ~ "Cue - No-Cue\n Δ cm/s",
        dv == "stride_mean" ~ "Δ cm/height",
        dv == "cadence_cv" ~ "Δ %",
        TRUE ~ dv
    )

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

    p <- ggplot(cue_vs_nocue, aes(x = subgroup, y = Est, color = contrast, group = contrast)) +
        geom_point(position = position_dodge(0.4), size = 4) +
        geom_errorbar(aes(ymin = CI_L, ymax = CI_U),
            position = position_dodge(0.4), width = 0.25
        ) +
        scale_x_discrete(labels = sample_sizes) +
        labs(
            title = paste0(title, " Interaction"),
            x = NULL,
            y = ylabel,
            color = "Cue Contrast"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.line = element_line(colour = "black")
        ) +
        scale_color_manual(values = c("#FF8C00", "#0072b2"))

    max_y <- max(cue_vs_nocue$CI_U, na.rm = TRUE)
    has_hc <- any(str_detect(m1c_int$subgroup_pairwise, "hc"))
    group_map <- if (has_hc) {
        c(hc = 1, nofog = 2, fog = 3)
    } else {
        c(nofog = 1, fog = 2)
    }

    sig_contrasts1 <- m1c_int %>%
        filter(
            str_detect(cue_pairwise, regex("no.?cue", ignore_case = TRUE)),
            str_detect(signif, "\\+")
        ) %>%
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
            y = (max_y) + max_y * .16 * seq_along(xmin)
        )

    if (nrow(sig_contrasts1) > 0) {
        ann_corrected <- sig_contrasts1 %>% filter(!is.na(signif.fdr) & signif.fdr != "")
        ann_uncorrected <- sig_contrasts1 %>% filter(is.na(signif.fdr) | signif.fdr == "")
        p <- p + geom_signif(
            data = ann_corrected,
            aes(xmin = xmin, xmax = xmax, annotations = signif.fdr, y_position = y, group = Est),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.04, textsize = 7.5, vjust = 0.5,
            color = "black", textcolor = "black"
        ) + geom_signif(
            data = ann_uncorrected,
            aes(xmin = xmin, xmax = xmax, annotations = signif, y_position = y, group = Est),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.04, textsize = 5.5, vjust = 0.1,
            color = "grey10"
        )
    }

    p
}

plots_int <- map(dv_list, ~ plot_int(df, dv = .x))
names(plots_int) <- dv_list

##################
### Save plots ###
##################
# set up boxplots
ymax <- 160
p_speed <- plots_box$speed_mean + coord_cartesian(xlim = c(0.5, 3.5), ylim = c(0, ymax), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, ymax, length.out = 5)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), ,
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

ymax <- 1.0
p_stride <- plots_box$stride_mean + coord_cartesian(xlim = c(0.5, 3.5), ylim = c(0, ymax), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, ymax, length.out = 5)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1) +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), ,
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

ymax <- 12
p_cadence <- plots_box$cadence_cv + coord_cartesian(xlim = c(0.5, 3.5), ylim = c(0, ymax), expand = FALSE) +
    scale_y_continuous(breaks = seq(0, ymax, length.out = 5)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1) +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), ,
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

# set up int plots
p2_speed <- plots_int$speed_mean + coord_cartesian(xlim = c(0.5, 3.5), ylim = c(-25, 25), expand = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1) +
    scale_y_continuous(breaks = seq(-25, 25, length.out = 5)) +
    theme(
        plot.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), ,
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

p2_stride <- plots_int$stride_mean + coord_cartesian(xlim = c(0.5, 3.5), ylim = c(-0.1, 0.1), expand = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1) +
    theme(
        plot.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), ,
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

p2_cadence <- plots_int$cadence_cv + coord_cartesian(xlim = c(0.5, 3.5), ylim = c(-3, 3), expand = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, linewidth = 1) +
    scale_y_continuous(breaks = seq(-3, 3, length.out = 5)) +
    theme(
        plot.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), ,
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )

# combine plots
c_plot <- (p_speed | p_stride | p_cadence) /
    (p2_speed | p2_stride | p2_cadence) +
    plot_layout(guides = "collect") &
    theme(
        text = element_text(size = 19),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 18)
    )

# Save plot
saveFile <- file.path(outputDir, "Figure-Kinematics_Contrast_6Panel.jpg")
ggsave(saveFile, plot = c_plot, width = 18, height = 9, dpi = 600)

print("done.")
