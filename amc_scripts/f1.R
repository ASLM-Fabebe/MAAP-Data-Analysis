#amc_prep
##adding folder_path

##load the packages
source('scripts/install_packages_packman.R')

if (rstudioapi::isAvailable()) {
  folder_path <- rstudioapi::selectDirectory()
  print(folder_path)
} else {
  #  cat("Not running in RStudio.\n")
  folder_path=NA
}

#importing files both xls or xlsx
input_file <- list.files(
  folder_path,
  pattern = "^AMC.*\\.xls[x]?$",
  ignore.case = TRUE
)

amc_raw <- readxl::read_excel(file.path(folder_path,input_file))

# Replace missing dates with dates e.g. from registration date col --------
excel_origin = "1899-12-30"
date_parse_vec = c("ymd HMS", "ymd HM", "mdy HMS", "dmy HMS", "ymd", "dmy", "mdy")


###matching the cols
choices1 <- c(names(amc_raw),'not available')  # Use your real variable


#prelim required columns
cols <- c("region", "product", "strength", "strength_unit",
          "strength_liquid_drugs", "volume_liquid_drugs",
          "pack_size","pack_size_unit", "quantity", "date", "route")

empty_amc_df <- data.frame(Required_variables=cols,
                           Corresponding_variables=c(rep('',11)))


# Match column classes of df2 to df1 and allow missing columns

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

#
bind_rows_match_classes <- function(dfs) {
  Reduce(function(x, y) bind_rows(x, match_col_classes(x, y)), dfs)
}

#create a folder to hold the temporary files
amc_updates_dir <- file.path(folder_path, "analysis_updates")

if(!dir.exists(amc_updates_dir)){dir.create(amc_updates_dir, recursive = T)}

#colors
my_colors <- c(
  brewer.pal(9, "Set1"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Dark2")
)

