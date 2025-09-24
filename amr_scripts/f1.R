
#load the packages
source('scripts/install_packages_packman.R')

##adding folder_path

if (rstudioapi::isAvailable()) {
  folder_path <- rstudioapi::selectDirectory()
  print(folder_path)
} else {
  #  cat("Not running in RStudio.\n")
  folder_path=NA
}


input_file <- list.files(folder_path, pattern = "^AMR.*.xlsx")

amr <- readxl::read_excel(file.path(folder_path,input_file), guess_max = 21474836) %>%
  dplyr::select(where(~ any(!is.na(.))))#drop all empty columns

# Replace missing dates with dates e.g. from registration date col --------
excel_origin = "1899-12-30"
date_parse_vec = c("ymd HMS", "ymd HM", "mdy HMS", "dmy HMS", "ymd", "dmy", "mdy")


###matching the cols when combining with reduce
# Sample vector
choices1 <- sort(c(names(amr),'not available'))  # Use your real variable

# Define initial df
initial_df <- data.frame(
  man_vars = c('Specimen date', 'Date of data entry', 'Specimen type', "Organism", 'Age','AST guidelines',
               'Sex', "Identification number", "Laboratory Name", "Patient Department", "Patient Location Type"),
  my_dataset = factor(rep('please select', 11), levels = choices1),
  stringsAsFactors = FALSE
)

#long_df parameters
long_df_cols <- data.frame(
  man_vars = c("Antibiotics_column", "AST_results","Identification number", "Organism", 'Specimen date'),
  my_dataset = factor(rep('please select', 5), levels = choices1),
   stringsAsFactors = FALSE
)


##Returning all columns present to allow user identify antibiotic columns
ab_cols_user <- tibble(my_dataset=unique(names(amr)), antibiotic_column=rep('', length(unique(names(amr)))))



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
#create the results directory
amr_updates_dir <- file.path(folder_path, "analysis_updates")

if(!dir.exists(amr_updates_dir)){dir.create(amr_updates_dir, recursive = T)}

#suppress summarize messages
options(dplyr.summarise.inform = FALSE)
