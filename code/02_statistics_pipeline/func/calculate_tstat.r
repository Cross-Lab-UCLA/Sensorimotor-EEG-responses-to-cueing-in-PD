# Function to calculate T-Stat in a dataframe with Est and SE.
library(dplyr)

calculate_tstat <- function(df, est_string, se_string) {
    if (!all(c(est_string, se_string) %in% names(df))) {
        stop(paste("Error: The input dataframe must contain columns named", est_string, "and", se_string))
    }

    df_calculated <- df %>%
        mutate(t_stat = !!sym(est_string) / !!sym(se_string))

    # insert t_stat column immediately after SE column
    se_index <- which(names(df_calculated) == se_string)
    all_cols <- names(df_calculated)
    pre_se_cols <- all_cols[1:(se_index)] # includes 'SE' and everything before it
    t_stat_col <- "t_stat"
    post_se_cols <- all_cols[(se_index + 1):length(all_cols)]
    final_column_order <- c(pre_se_cols, t_stat_col, post_se_cols)

    df_reordered <- df_calculated %>%
        select(all_of(final_column_order))

    # rounding
    df_final <- df_reordered %>%
        mutate(across(where(is.numeric), ~ round(.x, 4)))

    return(df_final)
}
