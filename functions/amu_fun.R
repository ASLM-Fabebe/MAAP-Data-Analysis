
# Match column classes of df2 to df1 and allow missing columns
match_col_classes <- function(df1, df2) {
  # Handle completely empty df2
  if (nrow(df2) == 0 && ncol(df2) == 0) {
    df2 <- as.data.frame(matrix(nrow = 0, ncol = 0))
  }

  # Add missing columns to df2
  missing_cols <- setdiff(names(df1), names(df2))
  for (col in missing_cols) {
    class_type <- class(df1[[col]])[1]
    df2[[col]] <- if (nrow(df2) == 0) {
      # special handling: zero-length vector of correct type
      if (class_type == "numeric") numeric(0) else
        if (class_type == "integer") integer(0) else
          if (class_type == "character") character(0) else
            if (class_type == "factor") factor(levels = levels(df1[[col]])) else
              logical(0)
    } else {
      # df2 has rows → fill with NA of correct type
      if (class_type == "numeric") rep(NA_real_, nrow(df2)) else
        if (class_type == "integer") rep(NA_integer_, nrow(df2)) else
          if (class_type == "character") rep(NA_character_, nrow(df2)) else
            if (class_type == "factor") factor(rep(NA, nrow(df2)), levels = levels(df1[[col]])) else
              rep(NA, nrow(df2))
    }
  }

  # Reorder columns to match df1
  df2 <- df2[, names(df1), drop = FALSE]

  # Match classes
  df2_matched <- map2_dfc(
    df2,
    df1,
    ~ {
      target_class <- class(.y)[1]
      if (target_class == "numeric") as.numeric(.x) else
        if (target_class == "integer") as.integer(.x) else
          if (target_class == "character") as.character(.x) else
            if (target_class == "factor") factor(.x, levels = levels(.y)) else
              .x
    }
  )

  names(df2_matched) <- names(df1)
  df2_matched
}



#Binding rows function
bind_rows_match_classes <- function(dfs) {
  Reduce(function(x, y) bind_rows(x, match_col_classes(x, y)), dfs)
}


##Dates
safe_as_posix <- function(x, ...) {
  out <- tryCatch(as.POSIXct(x, ...), error = function(e) NA)
  out
}


#processing date columns
##Dates
safe_as_posix <- function(x, ...) {
  tryCatch(as.POSIXct(x, ...), error = function(e) NA)
}

date_col_processing_vec <- function(x,
                                    excel_origin = as.Date("1899-12-30"),
                                    date_parse_vec = c("ymd", "mdy", "dmy")) {
  # raw input
  format_date_new <- x

  # numeric probe (for epochs + Excel serials)
  num <- suppressWarnings(as.numeric(format_date_new))

  # smart parsing to POSIXct
  posix <- dplyr::case_when(
    # already Date / POSIXt
    inherits(format_date_new, "Date")   ~ safe_as_posix(format_date_new, tz = "UTC"),
    inherits(format_date_new, "POSIXt") ~ safe_as_posix(format_date_new, tz = "UTC"),

    # numeric epochs: ms since 1970
    !is.na(num) & num > 1e12 ~ safe_as_posix(num / 1000,
                                             origin = "1970-01-01",
                                             tz = "UTC"),
    # numeric epochs: seconds since 1970
    !is.na(num) & num > 1e9  ~ safe_as_posix(num,
                                             origin = "1970-01-01",
                                             tz = "UTC"),

    # Excel serial days (anything numeric left over)
    !is.na(num) ~ safe_as_posix(num * 86400,
                                origin = excel_origin,
                                tz = "UTC"),

    # fallback: parse character dates in various formats
    TRUE ~ tryCatch(
      suppressWarnings(
        parse_date_time(as.character(format_date_new),
                        orders = date_parse_vec,
                        tz = "UTC")
      ),
      error = function(e) as.POSIXct(NA)
    )
  )

  # final cleaned Date
  as.Date(posix)
}




#colors for plotting
my_colors <- c(
  brewer.pal(9, "Set1"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Dark2")
)


##plotting frequencies of antibiotic uses

library(rlang)   # for sym()

plot_facility_prev_amu <- function(data, x1 = "ward_name", my_colors) {
  # turn the string into a column reference
  group_var <- sym(x1)

  facility_group_prev <- data %>%
    distinct(
      facility, !!group_var,
      date_of_data_collection, patient_id, antibiotic_names,
      .keep_all = TRUE
    ) %>%
    group_by(facility, year) %>%
    mutate(tot_facility = sum(n)) %>%
    ungroup() %>%
    group_by(facility, !!group_var, year) %>%
    summarise(
      tot_facility = mean(tot_facility),
      ward_count   = n(),
      .groups = "drop"
    ) %>%
    mutate(prev = round(ward_count * 100 / tot_facility, 1))

  write.csv(facility_group_prev, paste0(amu_dir,"/AMU_by_", x1,".csv"))


  if (nrow(facility_group_prev)>0) {

    p1 <- ggplot(
      facility_group_prev,
      aes(
        x = reorder(facility, -prev),
        y = prev,
        group = as.factor(year)
      )
    ) +
      geom_col(width = .9, aes(fill = !!group_var)) +
      labs(
        #title = paste("Percentage of AMU across", x1),
        x = "",
        y = "% of antibiotics used"
      ) +
      scale_fill_manual(values = my_colors) +
      theme_classic() +
      theme(
        legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 0)
      ) +
      coord_flip()
  }else{p1=ggplot(data.frame(x = c("a","b"), y = c(10, 20)), aes(x, y)) +
    geom_point() +
    labs(x = "", y = "") +
    annotate("text",
             x = "a", y = 15,
             label = "Not enough data to generate this plot") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  }

  ggsave(filename = paste0(amu_dir,"/AMU_by_",x1,".png"), plot = p1, width = 6, height = 6, dpi = 300)
}


