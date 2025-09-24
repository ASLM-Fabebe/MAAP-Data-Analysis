
#IF THE DATASET IS PROVIDED IN LONG FORMAT
cat('Checking if data reformatting is needed...\n')
message('Checking if data reformatting is needed...')

file_path <- file.path(amr_updates_dir,'long_format_amr_columns.xlsx')

if (file.exists(file_path)) {

  #read in the data
  long_cols <- read.xlsx(file_path)# or read.csv, fread, etc.
  long_cols$my_dataset[long_cols$man_vars=="Identification number"]
  #reshape the data
  amr <- amr %>% pivot_wider(id_cols = c(long_cols$my_dataset[long_cols$man_vars=="Identification number"],
                                  long_cols$my_dataset[long_cols$man_vars=="Organism"],everything()),
                      names_from = long_cols$my_dataset[long_cols$man_vars=="Antibiotics_column"],
                      values_from = long_cols$my_dataset[long_cols$man_vars=="AST_results"])

  cat("Data has been re-formatted..\n")
  message("Data has been re-formatted")
} else {
  amr <- amr

}


###call and rename stuff

cols_to_update <- read_excel(paste0(amr_updates_dir,"/select_amr_variables.xlsx"))

#move columns to beginning
# # Move given columns to the end of a data frame
# move_cols_to_beginning <- function(df, cols) {
#   all_cols <- names(df)
#   cols <- intersect(cols, all_cols)  # keep only valid ones
#   other_cols <- setdiff(all_cols, cols)
#   df <- df[, c(cols, other_cols), drop = FALSE]
#   return(df)
# }
#
# ##deduplicate names keep last
# dedup_cols_keep_last <- function(df) {
#   dupes <- duplicated(names(df), fromLast = TRUE)
#   df <- df[, !dupes, drop = FALSE]
#   return(df)
# }
#
# amr <- move_cols_to_beginning(amr, cols_to_update$man_vars)


# Perform renaming
names(amr)[names(amr) %in% cols_to_update$my_dataset] <- cols_to_update$man_vars[match(names(amr)[names(amr) %in% cols_to_update$my_dataset], cols_to_update$my_dataset)]

#Country code
#cntry=cols_to_update$enter_country_name_or_code[1]


#creating blank space holders for the unavailable cols
unavailable_cols=cols_to_update$man_vars[cols_to_update$my_dataset=='not available']

amr[unavailable_cols]=NA


get_inputs_and_prelim_cleanup <- function(){

  # Load AMR data sources ---------------------------------------------------

    safe_as_posix <- function(x, ...) {
     out <- tryCatch(as.POSIXct(x, ...), error = function(e) NA)
     out
   }

    amr <- amr %>%
     mutate(
       # pick source
       Specimen_date_new = if ("Date" %in% class(`Date of data entry`)) {
         coalesce(`Specimen date`, `Date of data entry`)
       } else {
         `Specimen date`
       },

       # numeric probe
       .num = suppressWarnings(as.numeric(Specimen_date_new)),

       # smart parsing
       .posix = dplyr::case_when(
         inherits(Specimen_date_new, "Date")    ~ safe_as_posix(Specimen_date_new, tz = "UTC"),
         inherits(Specimen_date_new, "POSIXt")  ~ safe_as_posix(Specimen_date_new, tz = "UTC"),

         # numeric epochs
         !is.na(.num) & .num > 1e12 ~ safe_as_posix(.num/1000, origin = "1970-01-01", tz = "UTC"), # ms
         !is.na(.num) & .num > 1e9  ~ safe_as_posix(.num,      origin = "1970-01-01", tz = "UTC"), # sec

         # Excel serial days
         !is.na(.num) ~ safe_as_posix(.num * 86400, origin = excel_origin, tz = "UTC"),

         # fallback to parsing strings
         TRUE ~ tryCatch(
           suppressWarnings(parse_date_time(as.character(Specimen_date_new),
                                            orders = date_parse_vec, tz = "UTC")),
           error = function(e) NA
         )
       ),

       specimen_date_cleaned = as.Date(.posix),

       # flag failures
       parse_failed = is.na(.posix) & !is.na(Specimen_date_new),
       parse_failed_value = ifelse(parse_failed, as.character(Specimen_date_new), NA_character_),

       r_id = row_number()
     ) %>%
     dplyr::select(-.num, -.posix)



   # Specify mandatory columns -----------------------------------------------

  man_cols = c("specimen_date_cleaned", # data specimen collected

               "Specimen type" , # type of specimemn

               "Organism"      # organism identified
               ) #date of collection

  # Remove columns with NA values in any of the mandatory columns

    #fill_function


  amr <- amr %>%

    filter(if_all(all_of(man_cols), ~ !is.na(.)))

  return(amr)

}



get_demographics <- function(df){
  # get demographics

  demo_vec <- c("r_id","Identification number","First name","Last name",

                "Sex","Date of birth","Age","Age category","Date of admission","Reason")

  lkp_demographics <- df %>%

    dplyr::select(any_of(demo_vec))

  return(lkp_demographics)
}



get_facilities_data <- function(df){

  # get lab and facility information

  facility_vec <- c("r_id","Identification number","Laboratory","Institution",

                    "Patient Location Type","Location type","Department","Origin","Country",

                    "Identification number", "Laboratory Name", "Patient Department")

  lkp_facility <- df %>%

    dplyr::select(any_of(facility_vec)) %>%

    #dplyr::select(where(~ !all(is.na(.))))  # remove columns with no values

  return(lkp_facility)
}


get_specimen_info <- function(df){
  # get specimens info

  specimen_vec <- c("r_id","Identification number","specimen_date_cleaned","Specimen number","Specimen type","Specimen type (Numeric)")

  lkp_specimens <-  df %>%

    dplyr::select(any_of(specimen_vec))

  return(lkp_specimens)

}

get_guideline_info <- function(df){
  # get specimens info

  ast_ref_vec <- c("r_id","AST guidelines")

  lkp_ast_ref <-  df %>%

    dplyr::select(any_of(ast_ref_vec))

  return(lkp_ast_ref)

}



file_path <- file.path(amr_updates_dir,'present_antibiotic_columns.xlsx')

if (file.exists(file_path)) {

  #read in the data
 ab_cols_user_vec <- read.xlsx(file_path) %>% # or read.csv, fread, etc.
  filter(antibiotic_column=='Yes')

 ab_cols_user_vec <- ab_cols_user_vec$my_dataset
} else {
  ab_cols_user_vec <- NULL

}



abx_vec_dict <- c(unique(readxl::read_excel(paste0('amr_resources/Antibiotic_Codes.xlsx'))$Code),
                  str_split_i(unique(readxl::read_excel(paste0('amr_resources/Antibiotic_Codes.xlsx'))$Code),'_',1),
                  unique(readxl::read_excel(paste0('amr_resources/Antibiotic_Codes.xlsx'))$AntiMicrobialAgent),
                  ab_cols_user_vec)


get_test_results <- function(df){

   # map antibiotics

  abx_conformed <- as.ab(abx_vec_dict)

  # Specify mandatory columns -----------------------------------------------
  man_cols = c("specimen_date_cleaned", # data specimen collected

               "Specimen type" , # type of specimemn

               "Organism",      # organism identified
               "specimen_date_cleaned") #date of collection

  abx_vec <- names(df)[tolower(names(df)) %in% tolower(c(abx_vec_dict))]


  # get antimicrobial results
  main_vars <- c("r_id","Identification number",man_cols)
  amr_res <- df %>%
    dplyr::select(any_of(c(main_vars,abx_vec)))


  #harmonizint the nomenclature
  p1 <- which(names(amr_res) == "Identification number")
  p2 <- which(names(amr_res) == "Organism")
  p3 <- which(names(amr_res) == "Specimen type")
  p4 <- which(names(amr_res) =="specimen_date_cleaned")

  names(amr_res)[p1] <- "uid"
  names(amr_res)[p2] <- "organism"
  names(amr_res)[p3] <- "specimen_type"
  names(amr_res)[p4] <- "specimen_date"



  amr_res <- amr_res %>%
    drop_na(organism, specimen_type,specimen_date )


  return(amr_res)

}


pivot_abx_results <- function(df){
  # Abx results preparation -------------------------------------------------

  abx_vec <- names(df)[tolower(names(df)) %in% tolower(c(abx_vec_dict))]

  amr_res <- df %>%

    mutate(bacteria = as.mo(organism, info = TRUE)) %>%

    mutate(gramstain = mo_gramstain(bacteria)) %>%

    dplyr::select(r_id, uid, specimen_type,bacteria, organism, gramstain,everything()) %>%

    mutate(across(any_of(abx_vec), ~as.character(.)))            #watch out for the uninterpretable bacteria and antibiotics

  #mo_failures()

  #mo_uncertainties()

  # pivot drugs from wide into long format and map to standard drug codes

  famr_long <- amr_res %>%

    tidyr::pivot_longer(names_to = "drug_code",

                        cols = any_of(abx_vec),

                        values_to = "vals") %>%

    mutate(ab = as.ab(drug_code))

  if(nrow(famr_long) > 0){

    famr_long=famr_long

  } else {
    famr_long=NULL
    #stop("famr_long seems empty.")
  }

  return(famr_long)

}


get_sir_interpr <- function(df){
  # Drugs with SIR interpretations already



  famr_long_sir <- df %>%left_join(lkp_ast_ref, by='r_id') %>%

    mutate(int_id=row_number()) %>%

    dplyr::filter(str_detect(vals,'(?i)R|I|S|SDD|NI')) %>%

    mutate(test_type=ifelse(grepl('_NM|_EM', drug_code), 'mic','disk'),

           # Use Gilbert's logic to determine if DISK or MIC

           guideline_inp=ifelse(grepl('_N', drug_code), 'CLSI',

                            ifelse(grepl('_E', drug_code), 'EUCAST', 'CLSI')),

           interpreted_res='',                    #to hold results for the next part
           guideline=`AST guidelines`,
           guideline= ifelse(is.na(guideline), guideline_inp, guideline),

           intrinsic_res_status=''

    )%>%
    filter(!is.na(bacteria))

  if(nrow(famr_long_sir) > 0){

    famr_long_sir=famr_long_sir

  } else {
    famr_long_sir=NULL
    #stop("famr_long seems empty.")
  }


  return(famr_long_sir)
}



get_con_interp <- function(df){
  # Drugs with breakpoints (MIC or DISK)

  famr_long_con <- df %>%left_join(lkp_ast_ref, by='r_id') %>%

    mutate(int_id=row_number()) %>%

    dplyr::filter(!str_detect(vals,'(?i)R|I|S|SDD|NI')) %>%

    # Use Gilbert's logic to determine if DISK or MIC

    mutate(test_type=ifelse(grepl('_NM|_EM', drug_code), 'mic','disk'),

           guideline_inp=ifelse(grepl('_N', drug_code), 'CLSI',

                                ifelse(grepl('_E', drug_code), 'EUCAST', 'CLSI')),

           interpreted_res='',                    #to hold results for the next part
           guideline=`AST guidelines`,
           guideline= ifelse(is.na(guideline), guideline_inp, guideline),

           intrinsic_res_status=''

    ) %>%
    filter(!is.na(bacteria))

  return(famr_long_con)
}

# Safe antibiogram wrapper
safe_antibiogram <- function(df, antimicrobials = NULL, min_n = 30) {
  tryCatch({
    if (!is.null(antimicrobials)) {
      out <- df %>% antibiogram(antimicrobials = antimicrobials, minimum = min_n)
    } else {
      out <- df %>% antibiogram(minimum = min_n)
    }

    # Return empty df if no rows
    if (is.null(out) || nrow(out) == 0) {
      return(data.frame())
    } else {
      return(out)
    }
  },
  error = function(e) {
    message("⚠ No valid antibiogram could be generated: ", e$message)
    return(data.frame())
  })
}


