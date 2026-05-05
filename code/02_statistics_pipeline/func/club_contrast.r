## club_contrast
# Function calls linear_contrast from ClubSandwich to get robust SE, CI, and p vals.
# It requires emmeans to get the contrast matrix for the contrasts. 
#
# Inputs:
#   - fitted model
#   - emmGrid from emmeans()
#       example: emmeans(model, ~ cue | subgroup)
#   - variance covariance matrix estimated using vcovCR
#       example: cov_cr <- vcovCR(model, type = "CR2", cluster = df$subject)
#
# LM 080725
##

## REQUIRED Libraries
library(clubSandwich)
library(emmeans)

## FUNCTION
club_contrast <- function(model, vcov_cr, emm, inter = FALSE) {
    if (inter) {
        contrasts <- contrast(emm, interaction = "pairwise")
    } else {
        contrasts <- pairs(emm)
    }

    contrast_mat <- contrasts@linfct
    grid_vars <- names(contrasts@grid) # grab all column names from contrasts

    results_list <- vector("list", nrow(contrast_mat))

    for (i in seq_len(nrow(contrast_mat))) {
        cArray <- contrast_mat[i, , drop = FALSE]

        res <- linear_contrast(
            obj = model,
            vcov = vcov_cr,
            contrasts = cArray,
            test = "Satterthwaite",
            p_values = TRUE
        ) %>% as.data.frame()

        ## get t-values for each contrast
        #res <- res %>%
        #mutate(t_statistic = Est / SE) %>% 
        #select(-p_val, p_val)   # move p_val to the end

        for (varname in grid_vars) {
            res[[varname]] <- contrasts@grid[[varname]][i]
        }

        results_list[[i]] <- res
    }

    out <- bind_rows(results_list)
    out <- out[, names(out) != "Coef"]
    out <- out[, c(grid_vars, setdiff(names(out), grid_vars))] %>%
        mutate(
            #signif = case_when(
            #    p_val < 0.001 ~ "***",
            #    p_val < 0.01 ~ "**",
            #    p_val < 0.05 ~ "*",
            #    p_val < 0.1 ~ ".",
            #    TRUE ~ ""
            #)
            signif = case_when(
                p_val < 0.05 ~ "+",
            )
        )
    rownames(out) <- NULL

    return(out)
}